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
struct GaussPoint {
    long double ksi;
    long double weight;
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
    double J1[2][2];
    double detJ;

    Jakobian() : detJ(0.0) {
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                J[i][j] = 0.0;
                J1[i][j] = 0.0;
            }
        }
    }

    void computeJacobian(const vector<Node>& elementNodes, const vector<long double>& dN_dEta, const vector<long double>& dN_dKsi) {
        J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;

        for (int i = 0; i < 4; ++i) {
            J[0][0] += dN_dEta[i] * elementNodes[i].x;
            J[0][1] += dN_dEta[i] * elementNodes[i].y;
            J[1][0] += dN_dKsi[i] * elementNodes[i].x;
            J[1][1] += dN_dKsi[i] * elementNodes[i].y;
        }

        detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

        if (detJ != 0) {
            J1[0][0] = J[1][1] / detJ;
            J1[0][1] = -J[0][1] / detJ;
            J1[1][0] = -J[1][0] / detJ;
            J1[1][1] = J[0][0] / detJ;
        }
    }

    void printJacobian() const {
        cout << fixed << setprecision(9);
        cout << "\nJacobian:" << endl;
        cout << "[ " << setw(12) << J[0][0] << " " << setw(12) << J[0][1] << " ]" << endl;
        cout << "[ " << setw(12) << J[1][0] << " " << setw(12) << J[1][1] << " ]" << endl;
        cout << "detJ= " << detJ << endl;
    }
};
struct IntegrationPointResults {
    Jakobian jakobian;
    vector<double> dN_dx;
    vector<double> dN_dy;
    vector<vector<double>> H_point;
};
struct Element {
    int ID[4];
    Jakobian jakobian;
    double Hbc[4][4];  // Existing boundary condition heat transfer matrix
    vector<double> localP;  // New vector to store local P for each element

    Element() {
        for (int i = 0; i < 4; ++i) {
            ID[i] = 0;
            for (int j = 0; j < 4; ++j) {
                Hbc[i][j] = 0.0;
            }
        }
        localP.resize(4, 0.0);  // Initialize localP with 4 elements set to 0
    }
};

struct GlobalData
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
    int npc; // liczba punktów całkowania
    vector<vector<long double>> dN_dEta; // rozmiar [npc][4]
    vector<vector<long double>> dN_dKsi; // rozmiar [npc][4]
    vector<double> N_i; // rozmiar [4]

    ElemUniv(int numPoints) : npc(numPoints), N_i(4, 0.0),
        dN_dEta(numPoints, vector<long double>(4, 0.0)),
        dN_dKsi(numPoints, vector<long double>(4, 0.0)) {
        // Konstruktor inicjalizuje wektory z odpowiednimi rozmiarami i wartościami zerowymi
    }

    void computeShapeFunctionDerivatives() {
        // Określenie liczby punktów Gaussa na jednym wymiarze (2 dla 2x2, 3 dla 3x3)
        int gaussOrder = static_cast<int>(sqrt(npc));

        // Sprawdzenie, czy npc to kwadrat liczby całkowitej
        if (gaussOrder * gaussOrder != npc) {
            throw std::invalid_argument("Invalid number of integration points (npc). Must be a perfect square (4, 9, etc.).");
        }

        // Generowanie współrzędnych punktów Gaussa
        vector<long double> gaussPoints(gaussOrder);
        if (gaussOrder == 2) {
            gaussPoints = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
        }
        else if (gaussOrder == 3) {
            gaussPoints = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
        }
        else {
            throw std::invalid_argument("Unsupported Gauss order. Supported orders are 2x2 and 3x3.");
        }

        // Generowanie wartości ksi i eta dla wszystkich punktów całkowania
        vector<long double> ksi_values(npc), eta_values(npc);
        for (int i = 0; i < gaussOrder; ++i) {
            for (int j = 0; j < gaussOrder; ++j) {
                ksi_values[i * gaussOrder + j] = gaussPoints[j];
                eta_values[i * gaussOrder + j] = gaussPoints[i];
            }
        }

        // Obliczanie pochodnych funkcji kształtu dla każdego punktu
        for (int i = 0; i < npc; ++i) {
            long double ksi = ksi_values[i];
            long double eta = eta_values[i];

            dN_dEta[i][0] = -0.25 * (1.0 - ksi);
            dN_dEta[i][1] = 0.25 * (1.0 - ksi);
            dN_dEta[i][2] = 0.25 * (1.0 + ksi);
            dN_dEta[i][3] = -0.25 * (1.0 + ksi);

            dN_dKsi[i][0] = -0.25 * (1.0 - eta);
            dN_dKsi[i][1] = -0.25 * (1.0 + eta);
            dN_dKsi[i][2] = 0.25 * (1.0 + eta);
            dN_dKsi[i][3] = 0.25 * (1.0 - eta);
        }
    }
    void printShapeFunctionDerivatives() const {
        cout << fixed << setprecision(9);

        cout << "        d N1/d Ksi      d N2/d Ksi      d N3/d Ksi      d N4/d Ksi" << endl;
        for (int i = 0; i < npc; ++i) {
            cout << "pc" << i + 1 << "  ";
            for (int j = 0; j < 4; ++j) {
                cout << setw(12) << dN_dEta[i][j] << " ";
            }
            cout << endl;
        }

        cout << "\n        d N1/d Eta     d N2/d Eta     d N3/d Eta     d N4/d Eta" << endl;
        for (int i = 0; i < npc; ++i) {
            cout << "pc" << i + 1 << "  ";
            for (int j = 0; j < 4; ++j) {
                cout << setw(12) << dN_dKsi[i][j] << " ";
            }
            cout << endl;
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
void readFile(string fileName, GlobalData& globalData, Grid** grid) {
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

void computeDerivativesDxDy(const Jakobian& jakobian,
    const vector<long double>& dN_dEta,
    const vector<long double>& dN_dKsi,
    vector<double>& dN_dx,
    vector<double>& dN_dy) {
    dN_dx.resize(4, 0.0);
    dN_dy.resize(4, 0.0);

    for (int i = 0; i < 4; ++i) {
        dN_dx[i] = jakobian.J1[0][0] * dN_dEta[i] + jakobian.J1[0][1] * dN_dKsi[i];
        dN_dy[i] = jakobian.J1[1][0] * dN_dEta[i] + jakobian.J1[1][1] * dN_dKsi[i];
    }
}


vector<IntegrationPointResults> computeIntegrationPointsResults(
    const vector<Node>& elementNodes,
    const ElemUniv& elemUniv,
    double conductivity,
    const vector<pair<GaussPoint, GaussPoint>>& points) {

    vector<IntegrationPointResults> results(elemUniv.npc);

   
    for (int i = 0; i < elemUniv.npc; i++) {
      
        results[i].jakobian.computeJacobian(elementNodes, elemUniv.dN_dEta[i], elemUniv.dN_dKsi[i]);

        
        computeDerivativesDxDy(results[i].jakobian,
            elemUniv.dN_dEta[i],
            elemUniv.dN_dKsi[i],
            results[i].dN_dx,
            results[i].dN_dy);

        
        results[i].H_point.resize(4, vector<double>(4, 0.0));
        for (int m = 0; m < 4; m++) {
            for (int n = 0; n < 4; n++) {
                results[i].H_point[m][n] = results[i].jakobian.detJ * conductivity *
                    (results[i].dN_dx[m] * results[i].dN_dx[n] +
                        results[i].dN_dy[m] * results[i].dN_dy[n]);
            }
        }
    }

    return results;
}


vector<vector<double>> computeFinalHMatrix(
    const vector<IntegrationPointResults>& results,
    const vector<pair<GaussPoint, GaussPoint>>& points) {

    vector<vector<double>> H_final(4, vector<double>(4, 0.0));

    for (int i = 0; i < results.size(); i++) {
        for (int m = 0; m < 4; m++) {
            for (int n = 0; n < 4; n++) {
                H_final[m][n] += results[i].H_point[m][n] *
                    points[i].first.weight *
                    points[i].second.weight;
            }
        }
    }

    return H_final;
}
void printDerivatives(const vector<IntegrationPointResults>& results) {
    cout << "\nPochodne dN/dx i dN/dy dla wszystkich punktow calkowania:" << endl;
    cout << "\td N1/d x\td N2/d x\td N3/d x\td N4/d x" << endl;
    for (int i = 0; i < results.size(); i++) {
        cout << "Pc" << i + 1 << ":\t";

        
        for (const auto& val : results[i].dN_dx) {
            cout << setw(12) << val << "\t";
        }
        cout << endl;

        
    }
    cout <<endl<< "\td N1/d y\td N2/d y\td N3/d y\td N4/d y" << endl;

    for (int i = 0; i < results.size(); i++) {
        cout << "Pc" << i + 1 << ":\t";


        for (const auto& val : results[i].dN_dy) {
            cout << setw(12) << val << "\t";
        }
        cout << endl;


    }
}
void printHMatrices(const vector<IntegrationPointResults>& results) {
    cout << "\nMacierze H dla poszczególnych punktow calkowania:" << endl;
    for (int i = 0; i < results.size(); i++) {
        cout << "\nPc " << i + 1 << ":" << endl;
        for (const auto& row : results[i].H_point) {
            for (const auto& val : row) {
                cout << setw(12) << val << " ";
            }
            cout << endl;
        }
    }
}

void printHMatrix(const vector<vector<double>>& H) {
    cout << "\nMatrix H:" << endl;
    for (const auto& row : H) {
        for (const auto& val : row) {
            cout << setw(12) << val << " ";
        }
        cout << endl;
    }
}
struct Surface {
    double pc[4][2]; // Gauss points coordinates 
    int npc; // Number of integration points 
    int nodeIds[2]; // Node IDs for this surface 
};
struct Solver {
    vector<vector<double>> globalH;
    vector<vector<double>> globalHbc; // New global Hbc vector
    vector<double> globalP; // New global P
    int matrixSize;

    Solver(int nodesNumber) : matrixSize(nodesNumber) {
        globalH.resize(matrixSize, vector<double>(matrixSize, 0.0));
        globalHbc.resize(matrixSize, vector<double>(matrixSize, 0.0));
        globalP.resize(matrixSize, 0.0);
    }

    void printGlobalH() const {
        cout << "\nGlobalna macierz H:" << endl;
        for (const auto& row : globalH) {
            for (const auto& val : row) {
                cout << setw(12) << val << " ";
            }
            cout << endl;
        }
    }

    void printGlobalHbc() const {
        cout << "\nGlobalna macierz Hbc:" << endl;
        for (const auto& row : globalHbc) {
            for (const auto& val : row) {
                cout << setw(12) << val << " ";
            }
            cout << endl;
        }
    }
    void printGlobalHSum() const {
        cout << "\nSuma globalnych macierzy H i Hbc:" << endl;
        for (int i = 0; i < matrixSize; ++i) {
            for (int j = 0; j < matrixSize; ++j) {
                cout << setw(12) << globalH[i][j] + globalHbc[i][j] << " ";
            }
            cout << endl;
        }
    }
    void printGlobalP() const {
        cout << "\nGlobalny wektor P:" << endl;
        for (const auto& val : globalP) {
            cout << setprecision(1) << val << " ";
        }
    }
};
void aggregateLocalHToGlobalH(vector<vector<double>>& globalH, const vector<vector<double>>& localH, const vector<int>& nodes) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            int globalI = nodes[i] - 1;
            int globalJ = nodes[j] - 1;
            globalH[globalI][globalJ] += localH[i][j];
        }
    }
}
void aggregateGlobalH(Grid* grid, const GlobalData& data, Solver& solver) {
    ElemUniv elemUniv(9);
    elemUniv.computeShapeFunctionDerivatives();

    
    vector<pair<GaussPoint, GaussPoint>> gaussPoints = gaussPoints2D(3);

    
    for (int elem = 0; elem < grid->ElementsNumber; elem++) {
        
        vector<Node> elementNodes;
        vector<int> nodeIndices;

        for (int i = 0; i < 4; i++) {
            int nodeIndex = grid->elements[elem].ID[i];
            nodeIndices.push_back(nodeIndex);
            elementNodes.push_back(grid->nodes[nodeIndex - 1]);  
        }

       
        vector<IntegrationPointResults> elemResults = computeIntegrationPointsResults(
            elementNodes, elemUniv, data.Conductivity, gaussPoints);

        
        vector<vector<double>> localH = computeFinalHMatrix(elemResults, gaussPoints);

        
        aggregateLocalHToGlobalH(solver.globalH, localH, nodeIndices);
    }
}
const double gaussPoints4[2][2] = {
    {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)},
    {1.0 / sqrt(3.0), -1.0 / sqrt(3.0)}
};

void calculateHbc(Element& element, const Grid* grid, const GlobalData& globalData, int gaussPointsCount) {
    // Reset Hbc matrix
    vector<vector<double>> Hbc(4, vector<double>(4, 0.0));

    // Define surfaces for the element
    Surface surfaces[4];

    // Define surfaces based on element node IDs
    // Surface points are in reference coordinates (ksi, eta)
    // Surface 1 (bottom): y = -1
    surfaces[0].npc = gaussPointsCount;
    surfaces[0].nodeIds[0] = element.ID[0];
    surfaces[0].nodeIds[1] = element.ID[1];
    for (int i = 0; i < gaussPointsCount; ++i) {
        surfaces[0].pc[i][0] = gaussPoints4[i][0];  // ksi varies
        surfaces[0].pc[i][1] = -1.0;                // eta is fixed at -1
    }

    // Surface 2 (right): x = 1
    surfaces[1].npc = gaussPointsCount;
    surfaces[1].nodeIds[0] = element.ID[1];
    surfaces[1].nodeIds[1] = element.ID[2];
    for (int i = 0; i < gaussPointsCount; ++i) {
        surfaces[1].pc[i][0] = 1.0;                 // ksi is fixed at 1
        surfaces[1].pc[i][1] = gaussPoints4[i][0];  // eta varies
    }

    // Surface 3 (top): y = 1
    surfaces[2].npc = gaussPointsCount;
    surfaces[2].nodeIds[0] = element.ID[2];
    surfaces[2].nodeIds[1] = element.ID[3];
    for (int i = 0; i < gaussPointsCount; ++i) {
        surfaces[2].pc[i][0] = -gaussPoints4[i][0]; // ksi varies (inverted for top surface)
        surfaces[2].pc[i][1] = 1.0;                 // eta is fixed at 1
    }

    // Surface 4 (left): x = -1
    surfaces[3].npc = gaussPointsCount;
    surfaces[3].nodeIds[0] = element.ID[3];
    surfaces[3].nodeIds[1] = element.ID[0];
    for (int i = 0; i < gaussPointsCount; ++i) {
        surfaces[3].pc[i][0] = -1.0;                // ksi is fixed at -1
        surfaces[3].pc[i][1] = -gaussPoints4[i][0]; // eta varies (inverted for left surface)
    }

    // Process each surface
    for (int surfaceIdx = 0; surfaceIdx < 4; ++surfaceIdx) {
        Surface& surface = surfaces[surfaceIdx];

        // Calculate Jacobian for surface integration
        double x1 = grid->nodes[surface.nodeIds[0] - 1].x;
        double y1 = grid->nodes[surface.nodeIds[0] - 1].y;
        double x2 = grid->nodes[surface.nodeIds[1] - 1].x;
        double y2 = grid->nodes[surface.nodeIds[1] - 1].y;
        double detJ = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2)) / 2.0;

        // Integrate Hbc for this surface
        for (int pc = 0; pc < surface.npc; ++pc) {
            double ksi = surface.pc[pc][0];
            double eta = surface.pc[pc][1];

            // Calculate shape functions
            double N1 = 0.25 * (1 - ksi) * (1 - eta);
            double N2 = 0.25 * (1 + ksi) * (1 - eta);
            double N3 = 0.25 * (1 + ksi) * (1 + eta);
            double N4 = 0.25 * (1 - ksi) * (1 + eta);

            // Shape function matrix
            double N[4] = { N1, N2, N3, N4 };

            // Compute Hbc contribution
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    Hbc[i][j] += globalData.Alfa *
                        N[i] * N[j] * detJ *
                        (pc < 2 ? 1.0 : 1.0); // Weight for 2 point Gauss quadrature
                }
            }
        }
    }

    // Copy Hbc to element
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            element.Hbc[i][j] = Hbc[i][j];
        }
    }
}

void aggregateHbc(Grid* grid, const GlobalData& data, Solver& solver) {
    for (int elem = 0; elem < grid->ElementsNumber; elem++) {
        // Compute Hbc for each element
        calculateHbc(grid->elements[elem], grid, data, 2);

        // Aggregate to global Hbc
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                int globalI = grid->elements[elem].ID[i] - 1;
                int globalJ = grid->elements[elem].ID[j] - 1;
                solver.globalHbc[globalI][globalJ] += grid->elements[elem].Hbc[i][j];
            }
        }
    }
}
void calculateLocalPVector(Element& element, const Grid* grid, const GlobalData& globalData, int gaussPointsCount) {
    // Reset local P vector
    element.localP.assign(4, 0.0);

    // Define surfaces for the element
    Surface surfaces[4];

    // Define surfaces based on element node IDs
    // Surface points are in reference coordinates (ksi, eta)
    // Surface 1 (bottom): y = -1
    surfaces[0].npc = gaussPointsCount;
    surfaces[0].nodeIds[0] = element.ID[0];
    surfaces[0].nodeIds[1] = element.ID[1];
    for (int i = 0; i < gaussPointsCount; ++i) {
        surfaces[0].pc[i][0] = gaussPoints4[i][0];  // ksi varies
        surfaces[0].pc[i][1] = -1.0;                // eta is fixed at -1
    }

    // Surface 2 (right): x = 1
    surfaces[1].npc = gaussPointsCount;
    surfaces[1].nodeIds[0] = element.ID[1];
    surfaces[1].nodeIds[1] = element.ID[2];
    for (int i = 0; i < gaussPointsCount; ++i) {
        surfaces[1].pc[i][0] = 1.0;                 // ksi is fixed at 1
        surfaces[1].pc[i][1] = gaussPoints4[i][0];  // eta varies
    }

    // Surface 3 (top): y = 1
    surfaces[2].npc = gaussPointsCount;
    surfaces[2].nodeIds[0] = element.ID[2];
    surfaces[2].nodeIds[1] = element.ID[3];
    for (int i = 0; i < gaussPointsCount; ++i) {
        surfaces[2].pc[i][0] = -gaussPoints4[i][0]; // ksi varies (inverted for top surface)
        surfaces[2].pc[i][1] = 1.0;                 // eta is fixed at 1
    }

    // Surface 4 (left): x = -1
    surfaces[3].npc = gaussPointsCount;
    surfaces[3].nodeIds[0] = element.ID[3];
    surfaces[3].nodeIds[1] = element.ID[0];
    for (int i = 0; i < gaussPointsCount; ++i) {
        surfaces[3].pc[i][0] = -1.0;                // ksi is fixed at -1
        surfaces[3].pc[i][1] = -gaussPoints4[i][0]; // eta varies (inverted for left surface)
    }

    // Process each surface
    for (int surfaceIdx = 0; surfaceIdx < 4; ++surfaceIdx) {
        Surface& surface = surfaces[surfaceIdx];

        // Calculate Jacobian for surface integration
        double x1 = grid->nodes[surface.nodeIds[0] - 1].x;
        double y1 = grid->nodes[surface.nodeIds[0] - 1].y;
        double x2 = grid->nodes[surface.nodeIds[1] - 1].x;
        double y2 = grid->nodes[surface.nodeIds[1] - 1].y;
        double detJ = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2)) / 2.0;

        // Integrate P for this surface
        for (int pc = 0; pc < surface.npc; ++pc) {
            double ksi = surface.pc[pc][0];
            double eta = surface.pc[pc][1];

            // Calculate shape functions
            double N1 = 0.25 * (1 - ksi) * (1 - eta);
            double N2 = 0.25 * (1 + ksi) * (1 - eta);
            double N3 = 0.25 * (1 + ksi) * (1 + eta);
            double N4 = 0.25 * (1 - ksi) * (1 + eta);

            // Shape function matrix
            double N[4] = { N1, N2, N3, N4 };

            // Compute local P contribution
            for (int i = 0; i < 4; ++i) {
                element.localP[i] += globalData.Alfa * globalData.Tot * N[i] * detJ *
                    (pc < 2 ? 1.0 : 1.0); // Weight for 2 point Gauss quadrature
            }
        }
    }
    for (int i = 0; i < 4; i++) {
        cout <<setprecision(1)<< element.localP[i] << " ";
    }
    cout << endl;
}

void aggregateLocalPToGlobalP(Grid* grid, Solver& solver, GlobalData globalData) {
    for (int elem = 0; elem < grid->ElementsNumber; elem++) {
        // Compute local P for each element
        calculateLocalPVector(grid->elements[elem], grid, globalData, 2);

        // Aggregate to global P
        for (int i = 0; i < 4; i++) {
            int globalI = grid->elements[elem].ID[i] - 1;
            solver.globalP[globalI] += grid->elements[elem].localP[i];
        }
    }
}




int main()
{
    GlobalData globalData;
    Grid* grid;
    string fileName = "test2_4_4.txt";
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
    elemUniv.printShapeFunctionDerivatives();
    
   
    vector<Node> nodes1;
    Node node;
    for (int i = 0; i < 4; i++)
    {
        nodes1.push_back(node);
    }
    double conductivity = 30;
    nodes1[0].x = 0.01;
    nodes1[0].y = -0.01;
    nodes1[1].x = 0.025;
    nodes1[1].y = 0.0;
    nodes1[2].x = 0.025;
    nodes1[2].y = 0.025;
    nodes1[3].x = 0.0;
    nodes1[3].y = 0.025;
    vector<pair<GaussPoint, GaussPoint>> points = gaussPoints2D(2);
    vector<IntegrationPointResults> results = computeIntegrationPointsResults(
        nodes1, elemUniv, conductivity, points);
    for (int i = 0; i < results.size(); i++) {
        cout << "\nPunkt całkowania " << i + 1 << ":" << endl;
        results[i].jakobian.printJacobian();
    }
    printDerivatives(results);
    printHMatrices(results);
    
    vector<vector<double>> H_final = computeFinalHMatrix(results, points);

    cout << "\nKońcowa macierz H:" << endl;
    for (const auto& row : H_final) {
        for (const auto& val : row) {
            cout << setw(12) << val << " ";
        }
        cout << endl;
    }
    cout << "\n\n=== OBLICZENIA DLA DANYCH Z PLIKU ===" << endl;

    
    ElemUniv elemUnivGrid(4);
    elemUnivGrid.computeShapeFunctionDerivatives();

    
    vector<pair<GaussPoint, GaussPoint>> gaussPoints = gaussPoints2D(2);

    
    for (int elem = 0; elem < grid->ElementsNumber; elem++) {
        cout << "\nElement " << elem + 1 << ":" << endl;

       
        vector<Node> elementNodes;
        for (int i = 0; i < 4; i++) {
            int nodeIndex = grid->elements[elem].ID[i] - 1; 
            elementNodes.push_back(grid->nodes[nodeIndex]);
        }

        
        vector<IntegrationPointResults> elemResults = computeIntegrationPointsResults(
            elementNodes, elemUnivGrid, globalData.Conductivity, gaussPoints);

        
        vector<vector<double>> elem_H_final = computeFinalHMatrix(elemResults, gaussPoints);

        cout << "\nMacierz H dla elementu " << elem + 1 << ":" << endl;
        printHMatrix(elem_H_final);
        
        
        /*for (int i = 0; i < elemResults.size(); i++) {
            cout << "\nPunkt całkowania " << i + 1 << " dla elementu " << elem + 1 << ":" << endl;
            elemResults[i].jakobian.printJacobian();
        }
        printDerivatives(elemResults);*/
    }
    Solver solver(globalData.NodesNumber);
    cout << "\n=== AGREGACJA MACIERZY H GLOBALNEJ ===" << endl;
    aggregateGlobalH(grid, globalData, solver);

    // Wyświetlenie globalnej macierzy H
    solver.printGlobalH();

    for (int i = 0; i < grid->ElementsNumber;i++) {
        calculateHbc(grid->elements[i], grid, globalData,2);

        // Print Hbc for each element
        cout << "Macierz Hbc dla elementu " << i+1 << ":\n";
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                cout << setw(10) << grid->elements[i].Hbc[j][k] << " ";
            }
            cout << endl;
        }
       
    }
    cout << "\n=== AGREGACJA MACIERZY HBC ===" << endl;
    aggregateHbc(grid, globalData, solver);

    // Print global Hbc
    solver.printGlobalHbc();
    solver.printGlobalHSum();
    cout << "\n=== OBLICZENIE I AGREGACJA WEKTORA P ===" << endl;
    aggregateLocalPToGlobalP(grid, solver, globalData );

    // Print global P vector
    solver.printGlobalP();
    return 0;
}
