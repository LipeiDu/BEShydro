//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include "../include/BEShydroLOGO.h"
#include <iostream>
#include <istream>
#include <fstream>

using namespace std;

void displayLogo(){
    
    cout << "=============================================================="   << endl;
    cout << "                                                              "   << endl;
    cout << "             _____ _____ _____ _         _                    "   << endl;
    cout << "            | __  |   __|   __| |_ _ _ _| |___ ___            "   << endl;
    cout << "            | __ -|   __|__   |   | | | . |  _| . |           "   << endl;
    cout << "            |_____|_____|_____|_|_|_  |___|_| |___|           "   << endl;
    cout << "                                  |___|                       "   << endl;
    cout << "                                                              "   << endl;
    cout << "  By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz   "   << endl;
    cout << "                                                              "   << endl;
    cout << "=============================================================="   << endl;
    
}

void displayCopyright(){
    
    cout << "                                                              "   << endl;
    cout << "  If you use the code, please cite:                           "   << endl;
    cout << "  (1) D. Bazow et al, arXiv: 1608.06577                       "   << endl;
    cout << "  (2) L. Du et al,                                            "   << endl;
    cout << "                                                              "   << endl;
    cout << "=============================================================="   << endl;
    
}
