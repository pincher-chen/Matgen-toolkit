# Space Group Transformation

## Process

### Get symmetric information

Process files to find information about space groups in files (`_symmetry_\name_H-M\name_Hall`), then convert the symmetric information into hall form and save it.

### Get transformation matrix

Processing symmetrical information to obtain two transformation matrices of coordinates. The transformation matrix **sym** holds symmetry about xyz, and the transformation matrix **offset** holds offset about xyz.


Examples:

* `x, y, z` $\Rightarrow$ $M_{sym}$: \[[1, 0, 0], [0, 1, 0], [0, 0, 1]]  $M_{offset}$: [0.0, 0.0, 0.0]

* `-x, 1/2+y, z` $\Rightarrow$  $M_{sym}$: \[[-1, 0, 0], [0, 1, 0], [0, 0, 1]]  $M_{offset}$: [0.0, 0.5, 0.0]


### Calculate coordinates

After the transformation matrix is obtained, the corresponding symmetry points can be calculated according to the transformation matrix.

For point A(x, y, z)
    $$A_{sym} = A * M_{sym} + M_{offset}$$

## Code
```
func(cif)
    // step 1 get symmetric information
    symm = cif.get_symm();
    name_HM = cif.get_name_HM();
    name_Hall = cif.get_name_Hall();

    if symm.empty == true
        if name_Hall.empty && name_HM.empty 
            symm.append({"x", "y", "z"});
        else if name_Hall.empty == false
            symm = SymOpsHall[name_Hall];
        else
            symm = SymOpsHall[HM2Hall[name_HM]];

    // step 2 get transformation matrix
    symm_arr = []
    offset_arr = []
    ex = ["-x", "-y", "-z"]
    for item in symm
        s = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        t = [0.0, 0.0, 0.0]

        for var i to 3 by i++ do
            if ex[i] in item[i]
                s[i][i] = -1.0
            else
                s[i][i] = 1.0

            t[i] = frac2double(item[i])
        
        symm_arr.append(s)
        offset_arr.append(t)
    end

    // step 3 calculate coordinates
    for cd in coordinates
        for var i to symm_arr.size by 
            new_cd = cd * symm_arr[i] + offset_arr[i]
end func
```