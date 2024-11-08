#include <iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include <iomanip>      
using namespace std;

struct Node {
    long double x, y;
    bool BC;

    Node() : x(0.0), y(0.0), BC(false) {}  // Default constructor
};

//struct Jakobian {
//    double J[2][2];
//    double J1[2][2];
//    double detJ;
//    vector<vector<double>> H;
//
//    Jakobian() {
//        for (int i = 0; i < 2; i++) {
//            detJ = 0.0;
//            for (int j = 0; j < 4; j++) {
//                J[i][j] = 0.0;
//                J1[i][j] = 0.0;
//            }
//        }
//    }
//    void computeJacobian(const vector<Node>& elementNodes, const long double dN_dEta[4][4], const long double dN_dKsi[4][4]) {
//        // Reset Jacobian matrix
//        J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
//
//        for (int i = 0; i < 4; ++i) {
//            J[0][0] += dN_dEta[0][i] * elementNodes[i].x;  
//            J[0][1] += dN_dEta[0][i] * elementNodes[i].y;
//            J[1][0] += dN_dKsi[0][i] * elementNodes[i].x;
//            J[1][1] += dN_dKsi[0][i] * elementNodes[i].y;
//        }
//
//        detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
//    }
//    void printJacobian(const long double dN_dEta[4][4], const long double dN_dKsi[4][4]) const {
//        cout << fixed << setprecision(9);
//
//        cout << "        d N1/d Ksi      d N2/d Ksi      d N3/d Ksi      d N4/d Ksi" << endl;
//        for (int i = 0; i < 4; ++i) {
//            cout << "pc" << i + 1 << "  ";
//            for (int j = 0; j < 4; ++j) {
//                cout << setw(12) << dN_dEta[i][j] << " ";
//            }
//            cout << endl;
//        }
//
//        cout << "\n        d N1/d Eta     d N2/d Eta     d N3/d Eta     d N4/d Eta" << endl;
//        for (int i = 0; i < 4; ++i) {
//            cout << "pc" << i + 1 << "  ";
//            for (int j = 0; j < 4; ++j) {
//                cout << setw(12) << dN_dKsi[i][j] << " ";
//            }
//            cout << endl;
//        }
//
//        
//            cout << "\nJacobian:"<< endl;
//            cout << "[ " << setw(12) << this->J[0][0] << " " << setw(12) << this->J[0][1] << " ]" << endl;
//            cout << "[ " << setw(12) << this->J[1][0] << " " << setw(12) << this->J[1][1] << " ]" << endl;
//            cout << "detJ= " << this->detJ << endl;
//       
//    }
//    
//};
struct Jakobian {
    double J[2][2];
    double J1[2][2]; // Odwrotność macierzy Jacobiana
    double detJ;
    vector<vector<double>> H;

    Jakobian() : detJ(0.0) {
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                J[i][j] = 0.0;
                J1[i][j] = 0.0;
            }
        }
        H.resize(4, vector<double>(4, 0.0)); // Inicjalizacja macierzy H
    }

    void computeJacobian(const vector<Node>& elementNodes, const long double dN_dEta[4][4], const long double dN_dKsi[4][4]) {
        // Resetowanie Jacobiana
        J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;

        for (int i = 0; i < 4; ++i) {
            J[0][0] += dN_dEta[0][i] * elementNodes[i].x;
            J[0][1] += dN_dEta[0][i] * elementNodes[i].y;
            J[1][0] += dN_dKsi[0][i] * elementNodes[i].x;
            J[1][1] += dN_dKsi[0][i] * elementNodes[i].y;
        }

        detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

        // Obliczanie odwrotności Jacobiana
        if (detJ != 0) {
            J1[0][0] = J[1][1] / detJ;
            J1[0][1] = -J[0][1] / detJ;
            J1[1][0] = -J[1][0] / detJ;
            J1[1][1] = J[0][0] / detJ;
        }
    }

    void computeDerivativesDxDy(const long double dN_dEta[4][4], const long double dN_dKsi[4][4], vector<vector<double>>& dN_dx, vector<vector<double>>& dN_dy) {
    dN_dx.resize(4, vector<double>(4, 0.0));
    dN_dy.resize(4, vector<double>(4, 0.0));

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            dN_dx[i][j] = J1[0][0] * dN_dEta[i][j] + J1[0][1] * dN_dKsi[i][j];
            dN_dy[i][j] = J1[1][0] * dN_dEta[i][j] + J1[1][1] * dN_dKsi[i][j];
        }
    }
}


    void computeHMatrix(const vector<vector<double>>& dN_dx, const vector<vector<double>>& dN_dy, double conductivity) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                H[i][j] = (dN_dx[0][i] * dN_dx[0][j] + dN_dy[0][i] * dN_dy[0][j]) * conductivity * detJ;
            }
        }
    }
    void printJacobian(const long double dN_dEta[4][4], const long double dN_dKsi[4][4]) const {
                cout << fixed << setprecision(9);
        
                cout << "        d N1/d Ksi      d N2/d Ksi      d N3/d Ksi      d N4/d Ksi" << endl;
                for (int i = 0; i < 4; ++i) {
                    cout << "pc" << i + 1 << "  ";
                    for (int j = 0; j < 4; ++j) {
                        cout << setw(12) << dN_dEta[i][j] << " ";
                    }
                    cout << endl;
                }
        
                cout << "\n        d N1/d Eta     d N2/d Eta     d N3/d Eta     d N4/d Eta" << endl;
                for (int i = 0; i < 4; ++i) {
                    cout << "pc" << i + 1 << "  ";
                    for (int j = 0; j < 4; ++j) {
                        cout << setw(12) << dN_dKsi[i][j] << " ";
                    }
                    cout << endl;
                }
        
                
                    cout << "\nJacobian:"<< endl;
                    cout << "[ " << setw(12) << this->J[0][0] << " " << setw(12) << this->J[0][1] << " ]" << endl;
                    cout << "[ " << setw(12) << this->J[1][0] << " " << setw(12) << this->J[1][1] << " ]" << endl;
                    cout << "detJ= " << this->detJ << endl;
               
            }
    void printDxDy(const vector<vector<double>>& dN_dx, const vector<vector<double>>& dN_dy) const {
        
                    cout << fixed << setprecision(9);
            
                    cout << "        d N1/d x      d N2/d x      d N3/d x      d N4/d x" << endl;
                    for (int i = 0; i < 4; ++i) {
                        cout << "pc" << i + 1 << "  ";
                        for (int j = 0; j < 4; ++j) {
                            cout << setw(12) << dN_dx[i][j] << " ";
                        }
                        cout << endl;
                    }
            
                    cout << "\n        d N1/d y     d N2/d y     d N3/d y     d N4/d y" << endl;
                    for (int i = 0; i < 4; ++i) {
                        cout << "pc" << i + 1 << "  ";
                        for (int j = 0; j < 4; ++j) {
                            cout << setw(12) << dN_dy[i][j] << " ";
                        }
                        cout << endl;
                    } 
    }

    void printHMatrix() const {
        cout << "\nMatrix H for each integration point:" << endl;
        for (const auto& row : H) {
            for (const auto& val : row) {
                cout << setw(12) << val << " ";
            }
            cout << endl;
        }
    }
};
struct Element {
    int ID[4];
    Jakobian jakobian;

    Element() {
        for (int i = 0; i < 4; ++i) ID[i] = 0;
    }
};

struct globalData
{
    int SimulationTime;
    int SimulationStepTime;
    int Conductivity;
    int Alfa;
    int Tot;
    int InitialTemp;
    int Density;
    int SpecificHeat;
    int NodesNumber;
    int ElementsNumber;
    int npc;
};
struct Grid {
    int NodesNumber;
    int ElementsNumber;
    Element* elements;
    Node* nodes;

    Grid() : NodesNumber(0), ElementsNumber(0), elements(nullptr), nodes(nullptr) {}

    void allocateMemory(int nodesNumber, int elementsNumber) {
        NodesNumber = nodesNumber;
        ElementsNumber = elementsNumber;
        elements = new Element[ElementsNumber];
        nodes = new Node[NodesNumber];
    }

    ~Grid() {
        delete[] elements;
        delete[] nodes;
    }
};

struct ElemUniv {
    int npc;
    long double dN_dEta[4][4];
    long double dN_dKsi[4][4];
    vector<double> N_i;

    ElemUniv(int numPoints) : npc(numPoints), N_i(4, 0.0) {
        for (int i = 0; i < npc; ++i) {
            for (int j = 0; j < 4; ++j) {
                dN_dEta[i][j] = 0.0;
                dN_dKsi[i][j] = 0.0;
            }
        }
    }

    void computeShapeFunctionDerivatives() {
        const long double ksi_values[4] = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0), -1.0 / sqrt(3.0) };
        const long double eta_values[4] = { -1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };

        for (int i = 0; i < npc; ++i) {
            long double ksi = ksi_values[i];
            long double eta = eta_values[i];

            dN_dEta[i][0] = -0.25 * (1.0 - eta);
            dN_dEta[i][1] = 0.25 * (1.0 - eta);
            dN_dEta[i][2] = 0.25 * (1.0 + eta);
            dN_dEta[i][3] = -0.25 * (1.0 + eta);

            dN_dKsi[i][0] = -0.25 * (1.0 - ksi);
            dN_dKsi[i][1] = -0.25 * (1.0 + ksi);
            dN_dKsi[i][2] = 0.25 * (1.0 + ksi);
            dN_dKsi[i][3] = 0.25 * (1.0 - ksi);
        }
    }
    

};
void printMatrix(const vector<vector<long double>>& matrix) {
    for (const auto& row : matrix) {
        for (const auto& value : row) {
            cout << setw(15) << value << " ";
        }
        cout << endl;
    }
}
void readFile(string fileName, globalData& globalData, Grid** grid) {
    ifstream file;
    file.open(fileName);
    if (!file.is_open())
    {
        cout << "Blad wczytywania: " << fileName << endl;
        return;
    }
    else {
        string line;
        while (!file.eof()) {
            getline(file, line);
            istringstream lineStream(line);
            string scan;
            lineStream >> scan;
            if (scan == "SimulationTime")
            {
                lineStream >> globalData.SimulationTime;
                cout << "SimulationTime " << globalData.SimulationTime << endl;
            }
            else if (scan == "SimulationStepTime")
            {
                lineStream >> globalData.SimulationStepTime;
                cout << "SimulationStepTime " << globalData.SimulationStepTime << endl;
            }
            else if (scan == "Conductivity")
            {
                lineStream >> globalData.Conductivity;
                cout << "Conductivity " << globalData.Conductivity << endl;
            }
            else if (scan == "Alfa")
            {
                lineStream >> globalData.Alfa;
                cout << "Alfa " << globalData.Alfa << endl;
            }
            else if (scan == "Tot")
            {
                lineStream >> globalData.Tot;
                cout << "Tot " << globalData.Tot << endl;
            }
            else if (scan == "InitialTemp")
            {
                lineStream >> globalData.InitialTemp;
                cout << "InitialTemp " << globalData.InitialTemp << endl;
            }
            else if (scan == "Density")
            {
                lineStream >> globalData.Density;
                cout << "Density " << globalData.Density << endl;
            }
            else if (scan == "SpecificHeat")
            {
                lineStream >> globalData.SpecificHeat;
                cout << "SpecificHeat " << globalData.SpecificHeat << endl;
            }
            else if (scan == "Nodesnumber")
            {
                lineStream >> globalData.NodesNumber;
                cout << "Nodes number " << globalData.NodesNumber << endl;
            }
            else if (scan == "Elementsnumber")
            {
                lineStream >> globalData.ElementsNumber;
                cout << "Elements number " << globalData.ElementsNumber << endl;
            }



            else if (scan == "*Node")
            {
                *grid = new Grid;
                (*grid)->allocateMemory(globalData.NodesNumber, globalData.ElementsNumber);
                int nodeID;
                char comma;
                string x, y;
                cout << "*Node" << endl;
                for (int i = 0; i < globalData.NodesNumber; i++)
                {
                    getline(file, line);
                    istringstream lineStream(line);
                    lineStream >> nodeID >> comma >> x >> y;
                    (*grid)->nodes[i].x = stold(x);
                    (*grid)->nodes[i].y = stold(y);
                    (*grid)->nodes[i].BC = false;
                    cout << "      " << nodeID << ". \tx=  " << setprecision(9) << (*grid)->nodes[i].x << ", \ty= " << setprecision(9) << (*grid)->nodes[i].y << endl;
                }
            }

            else if (scan == "*Element,")
            {
                cout << "*Element, type=DC2D4" << endl;
                for (int i = 0; i < globalData.ElementsNumber; i++)
                {
                    getline(file, line);
                    vector<int> elements;

                    istringstream lineStream(line);
                    string scan;

                    while (getline(lineStream, scan, ','))
                    {
                        scan.erase(scan.find_last_not_of(" \t") + 1);
                        scan.erase(0, scan.find_first_not_of(" \t"));

                        if (!scan.empty())
                        {
                            elements.push_back(stoi(scan));
                        }
                    }

                    for (int j = 1; j < 5; j++)
                    {
                        (*grid)->elements[i].ID[j - 1] = elements[j];
                    }

                    cout << "ID: " << i + 1 << "\t[ \t";

                    for (int j = 0; j < 4; j++)
                    {
                        cout << (*grid)->elements[i].ID[j];
                        if (j < 3) { cout << ", \t"; }
                    }
                    cout << "\t]" << endl;


                }
            }
            else if (scan == "*BC")
            {
                cout << "*BC" << endl;
                getline(file, line);
                vector<int> boundaryNodes;

                istringstream lineStream(line);
                string scan;
                while (getline(lineStream, scan, ','))
                {
                    scan.erase(scan.find_last_not_of(" \t") + 1);
                    scan.erase(0, scan.find_first_not_of(" \t"));

                    if (!scan.empty())
                    {
                        boundaryNodes.push_back(stoi(scan));
                    }
                }

                for (int i = 0; i < boundaryNodes.size(); i++)
                {
                    (*grid)->nodes[boundaryNodes[i] - 1].BC = true;
                    cout << " " << boundaryNodes[i] << ",  ";
                }
                cout << endl;
            }

        }


    }
}
struct GaussPoint {
    long double ksi;
    long double weight;
};

vector<GaussPoint> gaussPoints1D(int numPoints) {
    vector<GaussPoint> points;

    if (numPoints == 1) {
        points.push_back({ 0.0, 2.0 });  
    }
    else if (numPoints == 2) {
        points.push_back({ -1.0 / sqrt(3.0), 1.0 }); 
        points.push_back({ 1.0 / sqrt(3.0), 1.0 });
    }
    else if (numPoints == 3) {
        points.push_back({ -sqrt(3.0 / 5.0), 5.0 / 9.0 }); 
        points.push_back({ 0.0, 8.0 / 9.0 });
        points.push_back({ sqrt(3.0 / 5.0), 5.0 / 9.0 });
    }
    return points;
}

long double gaussIntegration1D(double (*f)(long double), int numPoints) {
    vector<GaussPoint> points = gaussPoints1D(numPoints);
    long double result = 0.0;

    for (const GaussPoint& p : points) {
        result += p.weight * f(p.ksi);  
    }
    return result;
}

vector<pair<GaussPoint, GaussPoint>> gaussPoints2D(int numPoints) {
    vector<GaussPoint> points1D = gaussPoints1D(numPoints);
    vector<pair<GaussPoint, GaussPoint>> points2D;

    for (const GaussPoint& p1 : points1D) {
        for (const GaussPoint& p2 : points1D) {
            points2D.push_back({ p1, p2 });
        }
    }
    return points2D;
}

long double gaussIntegration2D(double (*f)(long double, long double), int numPoints) {
    vector<pair<GaussPoint, GaussPoint>> points = gaussPoints2D(numPoints);
    long double result = 0.0;

    for (const auto& p : points) {
        result += p.first.weight * p.second.weight * f(p.first.ksi, p.second.ksi); 
    }
    return result;
}

double func1D(long double x) {
    return 5*x*x+3*x+6;  
}

double func2D(long double x, long double y) {
    return 5*x*x*y*y+3*x*y+6;  
}
int main()
{
    globalData globalData;
    Grid* grid;
    string fileName = "test1_4_4.txt";
    readFile(fileName, globalData, &grid);
    int numPoints = 2;
    cout << fixed << setprecision(25);
    long double result1D = gaussIntegration1D(func1D, 1);
    cout << "Calka w 1D wynosi 1 punktowy: "  << result1D << endl;
    long double result1D2 = gaussIntegration1D(func1D, 2);
    cout << "Calka w 1D wynosi 2 punktowy: " << result1D2 << endl;
    long double result1D3 = gaussIntegration1D(func1D, 3);
    cout << "Calka w 1D wynosi 3 punktowy: " << result1D3 << endl;
    long double result2D = gaussIntegration2D(func2D, 1);
    cout << "Calka w 2D wynosi 1 punktowy: " << result2D << endl;
    long double result2D2 = gaussIntegration2D(func2D, 2);
    cout << "Calka w 2D wynosi 2 punktowy: " << result2D2 << endl;
    long double result2D3 = gaussIntegration2D(func2D, 3);
    cout << "Calka w 2D wynosi 3 punktowy: " << result2D3 << endl;
    ElemUniv elemUniv(4);
    elemUniv.computeShapeFunctionDerivatives();
    vector<Node> nodes1; 
    Node node;
    for (int i = 0; i < 4; i++)
    {
        nodes1.push_back(node);
    }
    nodes1[0].x = 0.0;
    nodes1[0].y = 0.0;
    nodes1[1].x = 0.025;
    nodes1[1].y = 0.0;
    nodes1[2].x = 0.025;
    nodes1[2].y = 0.025;
    nodes1[3].x = 0.0;
    nodes1[3].y = 0.025;
    Jakobian jacobian;
    jacobian.computeJacobian(nodes1, elemUniv.dN_dEta, elemUniv.dN_dKsi);
    cout << "Element " << "1" << ":\n";
    jacobian.printJacobian(elemUniv.dN_dEta, elemUniv.dN_dKsi);
    jacobian.computeJacobian(nodes1, elemUniv.dN_dEta, elemUniv.dN_dKsi);
    cout << "Element " << "1:\n";
    jacobian.printJacobian(elemUniv.dN_dEta, elemUniv.dN_dKsi);

    vector<vector<double>> dN_dx, dN_dy;
    jacobian.computeDerivativesDxDy(elemUniv.dN_dEta, elemUniv.dN_dKsi, dN_dx, dN_dy);
    jacobian.printDxDy(dN_dx, dN_dy);

    double conductivity = 30;
    jacobian.computeHMatrix(dN_dx, dN_dy, conductivity);
    jacobian.printHMatrix();
        
    elemUniv.computeShapeFunctionDerivatives();

    
    
    return 0;
}
