#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <vector>

void run_tca(std::vector<std::vector<double>>& source, std::vector<int>& source_labels, 
             std::vector<std::vector<double>>& target, std::vector<int>& target_labels);

void run_bda(std::vector<std::vector<double>>& source, std::vector<int>& source_labels, 
             std::vector<std::vector<double>>& target, std::vector<int>& target_labels);

void run_coral(std::vector<std::vector<double>>& source, std::vector<int>& source_labels, 
               std::vector<std::vector<double>>& target, std::vector<int>& target_labels);


#endif