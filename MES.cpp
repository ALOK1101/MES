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
    vector<double> t; // wektor temperatur (rozwiązanie)
    int matrixSize;

    Solver(int nodesNumber) : matrixSize(nodesNumber) {
        globalH.resize(matrixSize, vector<double>(matrixSize, 0.0));
        globalHbc.resize(matrixSize, vector<double>(matrixSize, 0.0));
        globalP.resize(matrixSize, 0.0);
        t.resize(matrixSize, 0.0);
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
        }cout << endl;
    }
    vector<vector<double>> getAggregatedH() const {
        vector<vector<double>> aggregatedH(matrixSize, vector<double>(matrixSize, 0.0));
        for (int i = 0; i < matrixSize; ++i) {
            for (int j = 0; j < matrixSize; ++j) {
                aggregatedH[i][j] = globalH[i][j] + globalHbc[i][j];
            }
        }
        return aggregatedH;
    }
    void solveGauss() {
        vector<vector<double>> A = getAggregatedH(); // Macierz [H] + [Hbc]
        vector<double> b = globalP; // Wektor {P}

        int n = matrixSize;

        // Eliminacja współczynników
        for (int i = 0; i < n - 1; i++) {
            // Wybór elementu głównego
            int maxRow = i;
            double maxVal = abs(A[i][i]);

            for (int k = i + 1; k < n; k++) {
                if (abs(A[k][i]) > maxVal) {
                    maxVal = abs(A[k][i]);
                    maxRow = k;
                }
            }

            // Zamiana wierszy jeśli znaleziono lepszy element główny
            if (maxRow != i) {
                swap(A[i], A[maxRow]);
                swap(b[i], b[maxRow]);
            }

            // Eliminacja Gaussa
            for (int j = i + 1; j < n; j++) {
                double factor = A[j][i] / A[i][i];

                for (int k = i; k < n; k++) {
                    A[j][k] -= factor * A[i][k];
                }
                b[j] -= factor * b[i];
            }
        }

        // Rozwiązanie układu równań przez podstawienie wsteczne
        t[n - 1] = b[n - 1] / A[n - 1][n - 1];

        for (int i = n - 2; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += A[i][j] * t[j];
            }
            t[i] = (b[i] - sum) / A[i][i];
        }
    }

    // Metoda wyświetlająca wyniki (temperatury)
    void printResults() const {
        cout << "\nWyniki - temperatury w wezlach:" << endl;
        for (int i = 0; i < matrixSize; ++i) {
            cout << "Wezel " << i + 1 << ": " << fixed << setprecision(4) << t[i] << " C" << endl;
        }
    }

    // Metoda wyświetlająca układ równań
    void printEquationSystem() const {
        vector<vector<double>> H = getAggregatedH();
        cout << "\nUklad rownan [H + Hbc]{t} = {P}:" << endl;

        for (int i = 0; i < matrixSize; ++i) {
            cout << "[ ";
            for (int j = 0; j < matrixSize; ++j) {
                cout << setw(10) << fixed << setprecision(4) << H[i][j] << " ";
            }
            cout << "] [ t" << i + 1 << " ] = [ " << setw(10) << globalP[i] << " ]" << endl;
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

// Function to calculate shape functions for boundary conditions
vector<double> calculateShapeFunctions(double ksi) {
    vector<double> N(2);
    N[0] = (1.0 - ksi) / 2.0;
    N[1] = (1.0 + ksi) / 2.0;
    return N;
}

// Function to calculate Hbc for a surface
void calculateSurfaceHbc(double length, vector<vector<double>>& surfaceHbc, int numPoints, double alfa) {
    auto gaussPoints = gaussPoints1D(numPoints);

    for (const auto& point : gaussPoints) {
        vector<double> N = calculateShapeFunctions(point.ksi);

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                surfaceHbc[i][j] += alfa * N[i] * N[j] * length / 2.0 * point.weight;
            }
        }
    }
}

// Function to check if nodes form a boundary edge
bool isBoundaryEdge(const Node& node1, const Node& node2) {
    return node1.BC && node2.BC;
}

// Function to calculate surface length
double calculateSurfaceLength(const Node& node1, const Node& node2) {
    return sqrt(pow(node2.x - node1.x, 2) + pow(node2.y - node1.y, 2));
}

// Function to calculate local P vector for a surface
void calculateSurfaceP(double length, vector<double>& surfaceP, int numPoints, double alfa, double tot) {
    auto gaussPoints = gaussPoints1D(numPoints);

    for (const auto& point : gaussPoints) {
        vector<double> N = calculateShapeFunctions(point.ksi);

        for (int i = 0; i < 2; i++) {
            surfaceP[i] += alfa * tot * N[i] * length / 2.0 * point.weight;
        }
    }
}

// Function to calculate Hbc for an element
void calculateHbc(Element& element, Grid* grid, const GlobalData& data, int numPoints) {
    // Initialize Hbc matrix
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            element.Hbc[i][j] = 0.0;
        }
    }

    // Initialize local P vector
    element.localP = vector<double>(4, 0.0);

    // Define element edges (pairs of local node indices)
    vector<pair<int, int>> edges = { {0,1}, {1,2}, {2,3}, {3,0} };

    for (const auto& edge : edges) {
        int node1_idx = element.ID[edge.first] - 1;
        int node2_idx = element.ID[edge.second] - 1;

        Node& node1 = grid->nodes[node1_idx];
        Node& node2 = grid->nodes[node2_idx];

        if (isBoundaryEdge(node1, node2)) {
            // Calculate surface length
            double length = calculateSurfaceLength(node1, node2);

            // Calculate surface Hbc
            vector<vector<double>> surfaceHbc(2, vector<double>(2, 0.0));
            calculateSurfaceHbc(length, surfaceHbc, numPoints, data.Alfa);

            // Calculate surface P vector
            vector<double> surfaceP(2, 0.0);
            calculateSurfaceP(length, surfaceP, numPoints, data.Alfa, data.Tot);

            // Aggregate surface Hbc into element Hbc
            element.Hbc[edge.first][edge.first] += surfaceHbc[0][0];
            element.Hbc[edge.first][edge.second] += surfaceHbc[0][1];
            element.Hbc[edge.second][edge.first] += surfaceHbc[1][0];
            element.Hbc[edge.second][edge.second] += surfaceHbc[1][1];

            // Aggregate surface P into element localP
            element.localP[edge.first] += surfaceP[0];
            element.localP[edge.second] += surfaceP[1];
        }
    }
}

// Function to aggregate local Hbc matrices into global Hbc matrix
void aggregateHbc(Grid* grid, const GlobalData& data, Solver& solver) {
    for (int elem = 0; elem < grid->ElementsNumber; elem++) {
        Element& element = grid->elements[elem];

        // Aggregate Hbc
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                int globalI = element.ID[i] - 1;
                int globalJ = element.ID[j] - 1;
                solver.globalHbc[globalI][globalJ] += element.Hbc[i][j];
            }
        }
    }
}

// Function to aggregate local P vectors into global P vector
void aggregateLocalPToGlobalP(Grid* grid, Solver& solver, const GlobalData& data) {
    // Clear global P vector
    fill(solver.globalP.begin(), solver.globalP.end(), 0.0);

    // Aggregate local P vectors
    for (int elem = 0; elem < grid->ElementsNumber; elem++) {
        Element& element = grid->elements[elem];

        for (int i = 0; i < 4; i++) {
            int globalI = element.ID[i] - 1;
            solver.globalP[globalI] += element.localP[i];
        }
    }
}
void initializeTemperatures(Solver& solver, const GlobalData& data) {
    fill(solver.t.begin(), solver.t.end(), data.InitialTemp);
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

    cout << "\n=== OBLICZENIA HBC DLA ELEMENTÓW ===" << endl;
    for (int i = 0; i < grid->ElementsNumber; i++) {
        calculateHbc(grid->elements[i], grid, globalData, 2);

        cout << "\nMacierz Hbc dla elementu " << i + 1 << ":" << endl;
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                cout << setw(10) << grid->elements[i].Hbc[j][k] << " ";
            }
            cout << endl;
        }

        cout << "\nWektor P dla elementu " << i + 1 << ":" << endl;
        for (int j = 0; j < 4; j++) {
            cout << setw(10) << grid->elements[i].localP[j] << " ";
        }
        cout << endl;
    }

    cout << "\n=== AGREGACJA MACIERZY HBC ===" << endl;
    aggregateHbc(grid, globalData, solver);
    solver.printGlobalHbc();

    cout << "\n=== SUMA MACIERZY H + HBC ===" << endl;
    solver.printGlobalHSum();

    cout << "\n=== AGREGACJA WEKTORA P ===" << endl;
    aggregateLocalPToGlobalP(grid, solver, globalData);
    solver.printGlobalP();

    cout << "\n=== ROZWIAZYWANIE UKLADU ROWNAN ===" << endl;

    // Inicjalizacja temperatur początkowych
    initializeTemperatures(solver, globalData);

    // Wyświetlenie układu równań przed rozwiązaniem
    solver.printEquationSystem();

    // Rozwiązanie układu równań
    solver.solveGauss();

    // Wyświetlenie wyników
    solver.printResults();

    return 0;
}
