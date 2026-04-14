#include <iostream>
#include <vector>
#include "helper.h"
#include "algorithms.h"

using namespace std;

int main() {
    string source_file, target_file;
    int choice;

    cout << "========================================" << endl;
    cout << "  DOMAIN ADAPTATION TEST ENVIRONMENT" << endl;
    cout << "========================================" << endl;

    cout << "Enter Source Dataset filename (e.g., datasets/A.csv): ";
    cin >> source_file;
    cout << "Enter Target Dataset filename (e.g., datasets/Z.csv): ";
    cin >> target_file;

    vector<vector<double>> source, target;
    vector<int> source_labels, target_labels;

    cout << "\nLoading datasets..." << endl;
    load_dataset(source_file, source, source_labels);
    load_dataset(target_file, target, target_labels);
    cout << "Loaded Source (" << source.size() << " rows) and Target (" << target.size() << " rows).\n";

    while(true) {
        cout << "\n========================================" << endl;
        cout << "Select Algorithm:" << endl;
        cout << "1. TCA" << endl;
        cout << "2. BDA" << endl;
        cout << "3. CORAL" << endl;
        cout << "0. Exit" << endl;
        cout << "Choice: ";
        cin >> choice;

        if (choice == 0) break;
        else if (choice == 1) run_tca(source, source_labels, target, target_labels);
        else if (choice == 2) run_bda(source, source_labels, target, target_labels);
        else if (choice == 3) run_coral(source, source_labels, target, target_labels);
        else cout << "Invalid choice." << endl;
    }
    return 0;
}