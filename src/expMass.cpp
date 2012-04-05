/*
 * Copyright Krzysztof Miernik 2012
 * k.a.miernik@gmail.com 
 *
 * Distributed under GNU General Public Licence v3
 * Based on chartDrawer, prints table of Masses (or Sn, Sp, etc.)
 */

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <cmath>
#include "nuclide.h"

using namespace std;

// Translates Z into element name
//

Nuclide::Nuclide() {
    Z = 0;
    elementName = "none";
    N = 0;
    A = 0;
    massDefect = 0;
    massError = 0;
    halfLife = 0;
    halfLifeString = "none";
    extrapolated = false;
    spin = "none";
    primaryDecayMode = unknown;
}

void Nuclide::nameElement() {
    ifstream in("periodic.dat");
    if (in) {
        string line;
        int lineN = -1;
        while ( (lineN < Z)&&(getline(in,line)) ) {
                lineN++;
        }
        in.close();
        if (lineN < Z) {
            stringstream out;
            out << Z;
            line = "(" +  out.str() + ")";
        }
        elementName = line;
    } else {
            stringstream out;
            out << Z;
            elementName = "(" +  out.str() + ")";
    }
}

// Loads info from ame2003 entry (one line per isotope) into Nuclide struct
void process(Nuclide& t, const string& line) {
    /*
     * file nubtab03.asc is build using punch card scheme
     * certaint values are always found at certain positions in line thus
     * there's a lot of substrings below
     */
    // mass number
    t.A = atoi(line.substr(0,3).c_str());
    // atomic number and name of the element
    t.Z = atoi(line.substr(4,3).c_str());
    t.nameElement();
    // for convenience also neutron number is stored
    t.N = t.A - t.Z;
   
    // Mass defect
    string token = line.substr(18,9);
    // ...may include # if is based on extrapolation
    int found = token.find("#");
    if (found > -1) {
        token.erase(found);
        t.extrapolated = true;
    }
    else {
        t.extrapolated = false;
    }
    t.massDefect = atof(token.c_str());
    
    // Mass defect
    token = line.substr(29,9);
    // ...may include # if is based on extrapolation
    found = token.find("#");
    if (found > -1) 
        token.erase(found);
    t.massError = atof(token.c_str());

    // Half-life
    token = line.substr(60,7);
    // remove white spaces
    do {
        found = token.find(" ");
        if (found > -1)
            token.erase(found, 1);
    } while (found > -1); 
    // change plain '<' character into svg code 
    found = token.find("<");
    if (found > -1) {
        token.erase(found, 1);
        token.insert(found, "&lt; ");
    }
    // change plain '>' character into svg code 
    found = token.find(">");
    if (found > -1) {
        token.erase(found, 1);
        token.insert(found, "&gt; ");
    }
    // remove # characters (extrapolated values)
    found = token.find("#");
    if (found > -1) {
        token.erase(found, 1);
    }
    t.halfLifeString = token;
    // if half-life is existing then the next entry gives you units 
    // (ms, s, y, etc.)
    if (t.halfLifeString != "stbl" && t.halfLifeString != "p-unst") {
        token = line.substr(69,2);

        found = token.find(" ");
        if (found > -1)
            token.erase(found, 1);

        found = token.find(" ");
        if (found > -1)
            token.erase(found, 1);

        t.halfLifeString += " " + token;
    }
    // ground state spin value
    token = line.substr(79,13);
    // remove white spaces
    do {
        found = token.find(" ");
        if (found > -1)
            token.erase(found, 1);
    } while (found  > -1);
    // remove # 
    found = token.find("#");
    if (found > -1)
        token.erase(found);
    t.spin += token;
    
    // the rest of the lines gives decay modes
    // we are interested only in main one
    // we change all kinds of used characters (~, >, <, ?) into =
    token = line.substr(106);
    found = token.find("~");
    if (found > -1)
        token[found] = '=';
    found = token.find(">");
    if (found > -1)
        token[found] = '=';
    found = token.find("<");
    if (found > -1)
        token[found] = '=';
    found = token.find(" ?");
    if (found > -1)
        token[found] = '=';
    
    // Now the first entry is the largest branching
    found = token.find("=");
    if (found > -1) {
        token.erase(found);
    }
    else {
        found = token.find(";");
        if (found > -1)
            token.erase(found);
        do {
            found = token.find(" ");
            if (found > -1)
                token.erase(found,1);
        } while (found > -1);

        do {
            found = token.find("?");
            if (found > -1)
                token.erase(found,1);
        } while (found > -1);
    }

    // Store decay mode in struct
    if (t.halfLifeString == "stbl" || token == "IS")
        t.primaryDecayMode = stable;
    else if (t.halfLifeString == "p-unst")
        t.primaryDecayMode = unbound;
    else if (token == "B-")
        t.primaryDecayMode = betaM;
    else if (token == "B+" || token == "EC")
        t.primaryDecayMode = betaP;
    else if (token == "A")
        t.primaryDecayMode = alpha;
    else if (token == "SF")
        t.primaryDecayMode = fission;
    else if (token == "p")
        t.primaryDecayMode = proton;
    else if (token == "2p")
        t.primaryDecayMode = twoproton;
    else if (token == "n")
        t.primaryDecayMode = neutron;
    else if (token == "2n")
        t.primaryDecayMode = neutron;
    else
        t.primaryDecayMode = unknown;
}

int main() {
    ifstream nuFile("nubtab03km.asc");
    int sizeZ = 120; // Z = 0..119
    int sizeN = 178; // N = 0..177
    double **massArray= new double*[sizeZ];
    for (int i = 0; i < sizeZ; i++)
        massArray[i] = new double[sizeN]();

    double **errorArray= new double*[sizeZ];
    for (int i = 0; i < sizeZ; i++)
        errorArray[i] = new double[sizeN]();


    if (nuFile.good()) {
        string line;
        while (getline(nuFile, line) ) {
            if (line[0] != '#' && line[7] == '0') {
                Nuclide t;
                process(t, line);
                if (t.massDefect != 0) {
                //if (t.massDefect != 0 && t.extrapolated == false) {
                    massArray[t.Z][t.N] = t.massDefect / 1000.0;
                    errorArray[t.Z][t.N] = t.massError / 1000.0;
                }
                else if (t.massDefect == 0) {
                    massArray[t.Z][t.N] = 1e-12;
                    errorArray[t.Z][t.N] = t.massError / 1000.0;
               }

            }
        }
    }

    //cout << "# N  Z  S2p " << endl;
    cout << "# Z  N  S2n " << endl;

    // Change order of loops while changing S2p - S2n
    for (int i = 2; i < sizeZ; i++) {
    for (int j = 2; j < sizeN; j++) {
            //cout << i << " " << j << " " << massArray[i][j] << endl;
            // S2p
            /*
            if (massArray[i-2][j] != 0 && massArray[i][j] != 0)
                cout << j << " " << i << " " 
                     << 14.578 + massArray[i-2][j] - massArray[i][j]
                     << " " << sqrt(pow(errorArray[i-2][j],2) + pow(errorArray[i][j],2))
                     << endl;
            */
            // S2n
            if (massArray[i][j-2] != 0 && massArray[i][j] != 0)
                cout << i << " " << j << " " 
                     << 16.142 + massArray[i][j-2] - massArray[i][j]
                     << " " << sqrt(pow(errorArray[i][j-2],2) + pow(errorArray[i][j],2))
                     << endl;

        }
        cout << endl;
    }

    for (int i = 0; i < sizeZ; i++)
        delete []massArray[i];
    delete []massArray;

    return 0;
}
