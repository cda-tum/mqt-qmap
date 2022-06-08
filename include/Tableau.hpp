/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
*/

#ifndef QMAP_TABLEAU_HPP
#define QMAP_TABLEAU_HPP
#include <vector>


using tableau_row = std::vector<unsigned short>;

class Tableau {
public:
    Tableau() {
        tableau.clear();
    }

    std::vector<tableau_row> getTableau() {
        return tableau;
    }

    tableau_row operator[](std::size_t index){
        return tableau[index];
    }
private:
    std::vector<tableau_row> tableau{};
};


#endif //QMAP_TABLEAU_HPP
