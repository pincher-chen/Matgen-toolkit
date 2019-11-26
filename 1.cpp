#include <iostream> 

using namespace std;

void func(double lattice[][3]) {
    lattice[0][0] = -1;
    lattice[1][1] = 0;
    lattice[2][2] = 1;
}

struct S {
    double lattice[3][3];
};

int main() {
    S s;

    for(int i = 0; i < 3; i++) {
        for(int j  = 0; j < 3; j++) {
            s.lattice[i][j] = i+j;
        }
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cout << s.lattice[i][j] << " ";
        }
        cout << endl;
    }


    func(s.lattice);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cout << s.lattice[i][j] << " ";
        }
        cout << endl;
    }
}