#include "cctable.h"

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "\tUsage: " << argv[0] << " [vars file name]\n";
        
        return -1;
    }
    
    CcTable<unsigned int> input;
    
    if (!input.load_table(argv[1]))
    {
        return -1;
    }
    
    std::vector<int> positions;
    
    // Get all positions
    for (unsigned int i = 0; i < input.nrow(); i++)
    {
        bool found;
        std::string curr_pos = input.get_cell(i, 1, found);
        
        if (found)
        {
            positions.push_back(std::stoi(curr_pos));
        }
    }
    
    std::vector<int> positions_table;
    
    // Unique of the positions
    for (unsigned int p = 0; p < positions.size(); p++)
    {
        int curr_pos = positions.at(p);
        
        bool in_table = false;
        
        for (unsigned int t = 0; t < positions_table.size(); t++)
        {
            if (curr_pos == positions_table.at(t))
            {
                in_table = true;
                break;
            }
        }
        
        if (!in_table)
        {
            positions_table.push_back(curr_pos);
        }
    }
    
    // Count positions occurrencies
    std::vector<unsigned int> occurrencies;
    occurrencies.assign(positions_table.size(), 0);
    
    for (unsigned int o = 0; o < occurrencies.size(); o++)
    {
        unsigned int pos_to_check = positions_table.at(o);
        
        for (unsigned int p = 0; p < positions.size(); p++)
        {
            unsigned int curr_pos = (unsigned int)positions.at(p);
            
            if (curr_pos == pos_to_check)
            {
                occurrencies.at(o) += 1;
            }
        }
    }
    
    for (unsigned int s = 0; s < occurrencies.size(); s++)
    {
        std::cout << positions_table.at(s) << "\t" << occurrencies.at(s) << "\n";
    }
    
    return 0;
}