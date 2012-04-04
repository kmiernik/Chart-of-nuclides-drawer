/*
 * Krzysztof Miernik 22.08.2011
 * k.a.miernik@gmail.com
 */

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include "nuclide.h"

using namespace std;

Nuclide::Nuclide() {
    Z = 0;
    elementName = "none";
    N = 0;
    A = 0;
    massDefect = 0;
    halfLife = 0;
    halfLifeString = "none";
    extrapolated = false;
    spin = "none";
    primaryDecayMode = unknown;
}

// Translates Z into element name
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
    // we remove these addtional characters
    int found = token.find("#");
    if (found > -1)
        token.erase(found);
    t.massDefect = atof(token.c_str());
    
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
    ifstream nuFile("nubtab03.asc");
    // svg file header
    cout <<  "  <!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\"> " 
         << endl;
    cout << " <svg width=\"6200\" height=\"4000\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\"> " << endl;

    if (nuFile.good()) {
        string line;
        while (getline(nuFile, line) ) {
            if (line[0] != '#' && line[7] == '0') { //comment line & not isomer (which has [7] > 0)
                Nuclide t;
                process(t, line);
                unsigned x = t.N * 32;
                unsigned y = 4000 - (t.Z + 1) * 32;
                string rectStyle = "";
                string fontColor = "";
                switch (t.primaryDecayMode) {
                    case stable: rectStyle = "fill:#000000;"; fontColor = ";fill:#ffffff"; break;
                    case betaM: rectStyle = "fill:#758fff"; break;
                    case betaP: rectStyle = "fill:#ff7e75"; break;
                    case alpha: rectStyle = "fill:#fffe49"; break;
                    case fission: rectStyle = "fill:#5cbc57"; break;
                    case twoproton: 
                    case proton: rectStyle = "fill:#ffa425"; break;
                    case neutron: 
                    case unknown: 
                    case unbound: rectStyle = "fill:none;stroke-dasharray:2,2"; break;
                    default:rectStyle = "fill:none"; 
                }
                cout << " <rect style=\"stroke:#000000;stroke-width:0.5;" << rectStyle << "\" x=\""
                     << x << "\" y=\"" << y << "\" width=\"30\" height=\"30\"/>" << endl;
                cout << " <text";
                if (t.elementName[0] != '(') {
                    cout << " style=\"font-size:7px" << fontColor 
                         << "\" x=\"" << x + 12 - t.elementName.size() * 4;
                }
                else {
                    cout << " style=\"font-size:6px"
                         << "\" x=\"" << x + 12 - t.elementName.size() * 2 ;
                }
                cout << "\" y=\"" << y + 10 << "\">"  << endl
                     << t.elementName << t.A << endl
                     << "</text>" << endl;
                if (t.primaryDecayMode != stable && t.primaryDecayMode != unbound && t.primaryDecayMode != unknown) {
                    cout << " <text style=\"font-size:5px"
                    << fontColor << "\" x=\"" << x + 12 - t.halfLifeString.size()
                    << "\" y=\"" << y + 25 << "\">" << endl
                    << t.halfLifeString << endl << "</text>" << endl;
                cout << endl;
                }
            }
        }
    }
    cout << "</svg>" << endl;

    return 0;
}
