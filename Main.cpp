#include <iostream>
#include <vector>
#include "helper.h"
#include "algorithms.h"

using namespace std;

void load_datasets(string &source_file, string &target_file,
                   vector<vector<double>> &source, vector<vector<double>> &target,
                   vector<int> &source_labels, vector<int> &target_labels) {
    cout << "Enter Source Dataset filename (e.g., ../datasets/A.csv): ";
    cin >> source_file;
    cout << "Enter Target Dataset filename (e.g., ../datasets/Z.csv): ";
    cin >> target_file;

    source.clear(); target.clear();
    source_labels.clear(); target_labels.clear();

    load_dataset(source_file, source, source_labels);
    load_dataset(target_file, target, target_labels);
    cout << "Loaded Source (" << source.size() << " rows) and Target (" << target.size() << " rows).\n";
}







int main() {
    string source_file, target_file;
    int choice;

    cout << "======================================================" << endl;
    cout << "       DOMAIN ADAPTATION TEST ENVIRONMENT" << endl;
    cout << "======================================================" << endl;

    vector<vector<double>> source, target;
    vector<int> source_labels, target_labels;

    load_datasets(source_file, target_file, source, target, source_labels, target_labels);

    while (true) {
        cout << "\n========================================" << endl;
        cout << "Select Algorithm:" << endl;
        cout << "1. TCA" << endl;
        cout << "2. BDA" << endl;
        cout << "3. CORAL" << endl;
        cout << "4. Change Datasets" << endl;
        cout << "0. Exit" << endl;
        cout << "Choice: ";
        cin >> choice;

        if      (choice == 0) break;
        else if (choice == 1) run_tca(source, source_labels, target, target_labels);
        else if (choice == 2) run_bda(source, source_labels, target, target_labels);
        else if (choice == 3) run_coral(source, source_labels, target, target_labels);
        else if (choice == 4) load_datasets(source_file, target_file, source, target, source_labels, target_labels);
        else cout << "Invalid choice." << endl;
    }
    return 0;
}