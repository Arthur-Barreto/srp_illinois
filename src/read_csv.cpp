#include <fstream>
#include <iomanip> // For setting precision
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

struct Data {
    vector<double> time;
    vector<double> flux;
    vector<double> flux_err;
};

Data readCSV(const string &filename) {
    ifstream file(filename);
    Data data;
    string line;

    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return data;
    }

    // Read the header line
    if (getline(file, line)) {
        // Do nothing with the header line
    }

    while (getline(file, line)) {
        stringstream lineStream(line);
        string cell;
        vector<string> row;

        while (getline(lineStream, cell, ',')) {
            row.push_back(cell);
        }

        if (row.size() >= 3) { // Ensure there are at least 3 columns
            data.time.push_back(stod(row[0]));
            data.flux.push_back(stod(row[1]));
            data.flux_err.push_back(stod(row[2]));
        }
    }

    file.close();
    return data;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <filename>" << endl;
        return 1;
    }

    string filename = argv[1];
    Data data = readCSV(filename);

    // Set the decimal precision for output
    cout << fixed << setprecision(20);

    cout << "Time: ";
    for (const auto &value : data.time) {
        cout << value << " ";
    }
    cout << endl;

    cout << "Flux: ";
    for (const auto &value : data.flux) {
        cout << value << " ";
    }
    cout << endl;

    cout << "Flux_err: ";
    for (const auto &value : data.flux_err) {
        cout << value << " ";
    }
    cout << endl;

    return 0;
}
