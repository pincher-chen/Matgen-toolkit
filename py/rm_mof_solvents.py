#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import re
import math as m
import numpy as np
import collections
import argparse

parser = argparse.ArgumentParser(description='Handle MOFs cif and remove solvent atom(s)')
parser.add_argument('cif_in', help='input MOF cif file')
parser.add_argument('output_path', help='Output filepath')
parser.add_argument('-d',
                    '--skin_distance',
                    help='The skin distance(coefficient) you want to use (default=0.25)',
                    type=float,
                    default=0.25)
parser.add_argument('-s',
                    '--solvent',
                    help='Output the solvent was found',
                    action='store_true')
parser.add_argument('-m',
                    '--model',
                    help='Remove solvent molecules anyway',
                    action='store_true')

args = parser.parse_args()
filename = args.cif_in
output_path = args.output_path
distance = args.skin_distance
show = args.solvent
model = args.model
if distance < 1.0:
    skin_distance = distance
    coefficient = None
else:
    coefficient = distance
    skin_distance = None

class rmsolvent:
    def __init__(self,filename,outpath):
        self.filename = filename
        self.outpath = outpath

    def real_num(value):
        try:
            return float(value)
        except:
            if "(" in value:
                return float(value.split("(")[0])
            elif value in ['*','?']:
                return None
            else:
                print('return Value error, can not handle this!')
                exit(1)

    def readcif(self, filename):
        with open(filename,'r') as f:
            cf = f.read()
            #print(cf)
        return cf

    def split_cf(self, cf):
    ''' Split the information in the cif file and return as a dict '''
        total_loop_area = re.split('\s+?loop_\n', cf)
        loop_area = len(total_loop_area)
        loop_dict = {}
        for loop in total_loop_area:
            if '###########' in loop \
                    or re.match('data_.*?', loop)\
                    or re.match('\s*_symmetry_space_group\w+?', loop):
                chaos_loop = re.split('\n', loop)
                loop_dict['chaos_loop'] = chaos_loop
                #print(loop_dict)
            elif re.match('\s*_symmetry_equiv_pos\w+?', loop):
                sym_loop = re.split('\n', loop)
                loop_dict['sym_loop'] = sym_loop
                #print(sym_loop)
            elif re.match('\s*_atom_site\w+?', loop):
                site_loop = re.split('\n', loop)
                loop_dict['site_loop'] = site_loop
            elif re.match('\s*_atom_type\w+?', loop):
                type_loop = re.split('\n', loop)
                loop_dict['type_loop'] = type_loop
            elif re.match('\s*_geom_bond\w+?', loop):
                bond_loop = re.split('\n', loop)
                loop_dict['bond_loop'] = bond_loop
            else:
                unknow_loop = re.split('\n', loop)
                if unknow_loop is not None:
                    loop_dict['unknow_loop'] = unknow_loop
                else:
                    break
        #print(loop_dict)

        return loop_dict

    def site_loop(self, loop_dict):
    ''' Get the site information of the atom in the loop '''
        single_loop_dict = loop_dict.get('site_loop')
        #print(single_loop_dict)
        try:
            data_label, data_value = [], []
            for data_line in single_loop_dict:
                if re.match(r'\s*_\w+.*?', data_line):
                    new_data_line = data_line.replace(' ','')
                    data_label.append(new_data_line)
                else:
                    data_value.append(data_line)
        except:
            print('Failed to divide data')
            exit(1)
        else:
            index_range = len(data_label)
            loop_list = []
            for value_line in data_value:
                if len(value_line) != 0:
                    row = list(filter(None, value_line.split(' ')))
                else:
                    continue
                data_dict = {}
                for row_index in range(index_range):
                    data_key = data_label[row_index]
                    #print(data_key)
                    try:
                        row_value = row[row_index]
                    except IndexError:
                        continue
                    else:
                        if re.match(r'_atom_site_fract_\w',data_key):
                            data_dict[data_key] = rmsolvent.real_num(row_value)
                        else:
                            if '#END' == row_value:
                                break
                            else:
                                data_dict[data_key] = row_value
                            
                if len(data_dict) != 0:
                    loop_list.append(data_dict)
            #print(loop_list)

        return loop_list

    def seek_crystal(self, cf):
        ''' get crystal cell parameters from cif file'''
        cell_parameters_dict = {}
        cell_parameters_re = re.compile(r'(?P<name>_cell_[l|a].*?)\s+?(?P<value>.*)')
        cell_parameters = cell_parameters_re.findall(cf)
        for parameter in cell_parameters:
            cell_parameters_label, cell_parameters_value = parameter[0], \
                                                           rmsolvent.real_num(parameter[1])
            cell_parameters_dict[cell_parameters_label] = cell_parameters_value

        #print(cell_parameters_dict)
        return cell_parameters_dict

    def crystal(self, cell_parameters_dict):
        ''' calculate lattice contstant and angle'''
        try:
            a, b, c = cell_parameters_dict['_cell_length_a'], \
                          cell_parameters_dict['_cell_length_b'], \
                          cell_parameters_dict['_cell_length_c']
            ap, bt, ga = cell_parameters_dict['_cell_angle_alpha'] / 180 * m.pi, \
                             cell_parameters_dict['_cell_angle_beta'] / 180 * m.pi, \
                             cell_parameters_dict['_cell_angle_gamma'] / 180 * m.pi
        except KeyError:
            print('Cell parameter extract failed')
            exit(1)
        else:
            bc2 = b ** 2 + c ** 2 - 2 * b * c * m.cos(ap)

            h1 = a
            h2 = b * m.cos(ga)
            h3 = b * m.sin(ga)
            h4 = c * m.cos(bt)
            h5 = ((h2 - h4) ** 2 + h3 ** 2 + c ** 2 - h4 ** 2 - bc2) / (2 * h3)
            h6 = m.sqrt(c ** 2 - h4 ** 2 - h5 ** 2)
            self.cell = ([h1, 0., 0.], [h2, h3, 0.], [h4, h5, h6])
            self.angle = (ap, bt, ga)

        return self.cell, self.angle
        
    def get_atom_coordinates(self, site_loop_list):
    ''' get the coordinate information of the atom'''
        atom_cd_list = []
        for atom_dict in site_loop_list:
            #print(atom_dict)
            atom_cd_dict = {}
            atom_cd_x, atom_cd_y, atom_cd_z = atom_dict.get('_atom_site_fract_x'), \
                                              atom_dict.get('_atom_site_fract_y'), \
                                              atom_dict.get('_atom_site_fract_z')
            atom_cd = (atom_cd_x, atom_cd_y, atom_cd_z)
            species, atom_site = atom_dict.get('_atom_site_type_symbol'), \
                                 atom_dict.get('_atom_site_label')
            name = atom_site + '-' + species
            atom_cd_dict['label'], atom_cd_dict['atom_cd'] = name, atom_cd
            atom_cd_list.append(atom_cd_dict)
        #print(atom_cd_list)

        return atom_cd_list

    def get_cart(atom_array_cd, cell):
    '''Convert fractional coordinates to real coordinates and calculate bond lengths between different atoms'''
        lx_vector, ly_vector, lz_vector = cell[0], cell[1], cell[2]
        frac_a, frac_b, frac_c = atom_array_cd[0], \
                                 atom_array_cd[1], \
                                 atom_array_cd[2]
        x_cart = lx_vector[0] * frac_a + ly_vector[0] * frac_b + lz_vector[0] * frac_c
        y_cart = lx_vector[1] * frac_a + ly_vector[1] * frac_b + lz_vector[1] * frac_c
        z_cart = lx_vector[2] * frac_a + ly_vector[2] * frac_b + lz_vector[2] * frac_c

        return [x_cart, y_cart, z_cart]

    def calc_distance(a_array, b_array, cell):
        a_cart, b_cart = rmsolvent.get_cart(a_array, cell), \
                         rmsolvent.get_cart(b_array, cell)
        delt_d = np.array(a_cart) - np.array(b_cart)

        distance = m.sqrt(delt_d[0] ** 2 + delt_d[1] ** 2 + delt_d[2] ** 2)

        return distance

    def bond_distance(self, atom_cd_list, cell, angel):
        tt_atom_num = len(atom_cd_list)
        #print(tt_atom_num)
        dt_list = []
        for i in range(tt_atom_num):
            rest_atom_pos_list = atom_cd_list.copy()
            rest_atom_pos_list.remove(atom_cd_list[i])
            for atom_b in rest_atom_pos_list:
                dt_dict = {}
                atom_a = atom_cd_list[i]
                a, b = np.array(atom_a.get('atom_cd')),\
                        np.array( atom_b.get('atom_cd'))
                distance = rmsolvent.calc_distance(a, b, cell)
                a_label, b_label = atom_a.get('label'), atom_b.get('label')
                dt_name = a_label + '-' + b_label
                dt_dict[dt_name] = distance
                dt_list.append(dt_dict)
        #print(dt_list)
        return dt_list

    @staticmethod
    def atom_radius():
        radius_list = []
        with open('./conf/raduis.txt', 'r') as f:
            data = f.readlines()
        for line in data:
            radius_dict = {}
            radius = re.split('\s+?', line)
            chem_species, data = radius[0], radius[1:-1]
            radius_dict['species'] = chem_species
            radius_dict['radius'] = data[1]
            radius_dict['O_metal'] = data[-1]
            radius_list.append(radius_dict)
        # print(radius_list)

        return radius_list

    @staticmethod
    def known_solvent():
        solvent_list = []
        with open('./conf/solvent.txt', 'r') as f:
            data = f.readlines()
        for line in data:
            solvent_dict = {}
            solvent = re.split('\s+', line)
            # attention
            chem_formula, elements = solvent[0], solvent[1:-1]
            solvent_dict['chem_formula'] = chem_formula
            solvent_dict['elements'] = elements
            solvent_list.append(solvent_dict)
         # print(solvent_list)
        return solvent_list

    def judge_if_have_bond(self, dt_list, solvent_list, radius_list, skin_distance):
        bond = False
        have_bond_list = []
        for dt in dt_list:
            dt_mof = rmsolvent.get_cif_bond(dt)
            #print(dt_mof)
            atom_a, atom_b, dt_value = dt_mof[0], dt_mof[1], dt_mof[2]
            for radius_dict in radius_list:
                atom_species = rmsolvent.get_atom_radius(radius_dict)[0]
                atom_radius = rmsolvent.get_atom_radius(radius_dict)[1]
                if atom_a == atom_species: radius_a = float(atom_radius)
                if atom_b == atom_species: radius_b = float(atom_radius)
            if skin_distance is not None:
                bond_ab = dt_value - radius_a - radius_b
                if bond_ab < skin_distance:
                    bond = True
                    have_bond_list.append(dt)
            else:
                #print(radius_a)
                bond_ab = (radius_a + radius_b) * coefficient
                if dt_value < bond_ab:
                    bond = True
                    have_bond_list.append(dt)
                #print('have bond ', dt)
        #print(have_bond_list)
        print('The number of bonded atom pairs is ',len(have_bond_list))

        return bond, have_bond_list

    def get_cif_bond(dt_dict):
        for bond_type, value in dt_dict.items():
            atom_a, atom_b = bond_type.split('-')[1],  bond_type.split('-')[3]
            distance = value

        return atom_a, atom_b, distance 

    def get_atom_radius(radius_dict):
        atom = radius_dict.get('species')
        radius = radius_dict.get('radius')

        return atom, radius

    def same_atom(have_bond_dict):
        for bond_type, value in have_bond_dict.items():
            atom_a_site, atom_b_site = bond_type.split('-')[0],  bond_type.split('-')[2]

        return atom_a_site, atom_b_site

    def atom_subset(connect_list):
        ''' Grouping bonded atoms '''
        cn_atom_chain = []
        pair_num = len(connect_list)
        for i in range(pair_num):
            for j in range(pair_num):
                x = list(set(connect_list[i] + connect_list[j]))
                y = len(connect_list[j])+len(connect_list[i])
                if i == j or connect_list[i] == 0 or connect_list[j] == 0:
                    break
                elif len(x) < y:
                    connect_list[i] = x
                    connect_list[j] = [0]
        #print(connect_list)
        for chain in connect_list:
            if chain != [0]:
                cn_atom_chain.append(chain)
        #print(cn_atom_chain)
        return cn_atom_chain

    def connect_network(self, have_bond_list, atom_cd_fra):
        #print(have_bond_list)
        connect_list, cn_species = [], []
        for bond_dict in have_bond_list:
            atom_a_site, atom_b_site = rmsolvent.same_atom(bond_dict)[0],\
                                       rmsolvent.same_atom(bond_dict)[1]
            atom_pair = [atom_a_site,atom_b_site]
            connect_list.append(atom_pair)
        cn_atom_chain = rmsolvent.atom_subset(connect_list)
        #print(cn_atom_chain)
        for atom_chain in cn_atom_chain:
            cn_chain_dict = {}
            #print(atom_chain)
            per_cn_species = rmsolvent.get_cn_atom_species(atom_chain, atom_cd_fra)
            per_cn_species.sort()
            cn_chain_dict['atom_site'], cn_chain_dict['species']  = atom_chain, per_cn_species
            #print(cn_chain_dict)
            cn_species.append(cn_chain_dict)
        
        #print(cn_species)
        return cn_species

    def get_cn_atom_species( atom_chain, atom_cd_fra):
        ''' Get the atomic element symbol '''
        per_cn_species = []
        #print(atom_cd_fra)
        for atom_dict in atom_cd_fra:
            for atom_site_label in atom_chain:
                site_label, type_symbol = atom_dict['label'].split('-')[0],\
                                          atom_dict['label'].split('-')[1]
                if atom_site_label == site_label:
                    atom_species = type_symbol
                    per_cn_species.append(atom_species)
                    #print(per_cn_species)
                else:
                    continue
        #print(per_cn_species)
        return per_cn_species

    def cmp_solvent_known(self, cn_species, solvent_list, radius_list):
        ''' Compare existing list of solvents'''
        chem_list = rmsolvent.get_chem_formula(cn_species, radius_list) 
        print('The calculated solvent molecule to be screened is', chem_list[1])
        print('The MOF framework is ', chem_list[0])
        solvent_formula_list = []
        for solvent in solvent_list:
            solvent_formula = solvent.get('chem_formula')
            solvent_formula_list.append(solvent_formula)
        #print(solvent_formula_list)
        chk_solvent = chem_list[1]
        inlist = []
        for unknow_chem in chk_solvent:
            if unknow_chem in solvent_formula_list:
                print('Found ' +  unknow_chem +\
                      '\nStarting to remove...')
                inlist.append('1')
            else:
                print('Unable to find information about this molecule '\
                      +  unknow_chem)
                inlist.append('0')

        return chk_solvent,inlist 
 
    def cn_counter(chain):
        ''' Get the chemical formula of the atoms in the group'''
        chain_species = chain.get('species')
        counter = collections.Counter(chain_species)
        chem_formula = ''
        for chem,index in counter.items():
            atom_state = radius_list
            ele_n = chem + str(index)
            chem_formula += ele_n

        return chem_formula

    def get_chem_formula(cn_species, radius_list):
        framwork_list, solvent_chk_list = [], []
        for chain in cn_species:
            chem_formula = rmsolvent.cn_counter(chain)
            state_list = []
            chain_species = chain.get('species')
            for element in chain_species:
                state =  rmsolvent.get_atom_state(element, radius_list)
                state_list.append(state)
            if '1' in state_list:
                framwork_list.append(chem_formula)
            else:
                solvent_chk_list.append(chem_formula)

        return framwork_list, solvent_chk_list

    def get_atom_state(element, radius_list):
        '''return the element is metal or not'''
        for atom in radius_list:
            atom_species = atom.get('species')
            if element == atom_species:
                state = atom.get('O_metal')
                #print(state)
        if state is not None:
            return state
        else:
            print('Unable to find the state information of the corresponding element...')
            exit(1)

    def modify_cif_file(self, cf, chk_solvent, cn_species, output_path, filename):
        rm_atom_site = []
        inlist = chk_solvent[1]
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        for chain in cn_species:
            chem_formula = rmsolvent.cn_counter(chain)
            if chem_formula in chk_solvent[0]:
                atom_site = chain.get('atom_site')
                rm_atom_site.append(atom_site)
            else:
                continue
        #print(rm_atom_site)
        cf_file = os.path.split(filename)[1]
        #print(inlist)
        if len(rm_atom_site) != 0:
            if '0' not in inlist:
                clean_name = str(cf_file).split('.')[0] + '_chk.cif' 
            else:
                clean_name = str(cf_file).split('.')[0] + '_ndchk.cif'
            with open(output_path + os.sep + clean_name, 'w') as f:
                rmlist_index = []
                list_range = len(inlist)
                for num in range(list_range):
                    if inlist[num] == '0':
                        rmlist_index.append(num)
                #print(rmlist_index)
                if model:
                    print('Comparing known solvent lists, do not delete if this molecule don\'t exist\
                          \nCreating *ndchk.cif')
                    chk_list = []
                    for index in rmlist_index:
                        item = rm_atom_site[index]
                        chk_list.append(item)
                    for chk_item in chk_list:
                        rm_atom_site.remove(chk_item)
                else:
                    print('Directly delete calculated solvent molecules\
                           \nCreating *ndchk.cif')
                cf_list = re.split('\n', cf)
                #print(len(cf_list))
                for single_solvent in rm_atom_site:
                    del_line_list = []
                    for atomsite in single_solvent:
                        for line in cf_list:
                            if re.match('^' + atomsite + '\s+?', line):
                                #print(line)
                                cf_list.remove(line)
                            else:
                                continue
                #print(len(cf_list))
                for line in cf_list:
                    f.writelines(line + '\n')
        else:
            print('Solvent molecules not found, no modification required...')
            exit(0)                

if __name__ == '__main__':
    print('\nLooking for solvent in ', filename)
    go = rmsolvent(filename,output_path)
    cf = go.readcif(filename)
    loop_dict = go.split_cf(cf)
    cell_parameters_dict = go.seek_crystal(cf)
    cell_parameter = go.crystal(cell_parameters_dict)
    #print(cell_parameter)
    cell,angle = cell_parameter[0], cell_parameter[1]
    site_loop = go.site_loop(loop_dict)
    #print(site_loop)
    solvent_list, radius_list = go.known_solvent(),go.atom_radius()
    site_list = go.get_atom_coordinates(site_loop)
    dt_list = go.bond_distance(site_list, cell, angle)
    bond = go.judge_if_have_bond(dt_list,solvent_list,radius_list,skin_distance)
    have_bond_list = bond[1]
    cn_species = go.connect_network(have_bond_list, site_list)
    chk_solvent = go.cmp_solvent_known(cn_species,solvent_list, radius_list)
    if not show:
        go.modify_cif_file(cf,chk_solvent,cn_species,output_path,filename)
    else:
        exit(0)
