#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <cstdio>

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cerr << "Usage:\t" << argv[0] << " [file] [lines to be picked]\n";
        return 1;
    }
    
    srand(time(0));
    std::ifstream input;
    input.open(argv[1]);
    
    unsigned int lines = atoi(argv[2]);
    
    if (!input.is_open())
    {
        std::cerr << "Unable to open " << argv[0] << "\n";
        
        return 2;
    }
    
    char curr_line[8192];
    
    std::vector<std::string> all_lines;
    
    while (input.good())
    {
        input.getline(curr_line, 8192);
        
        std::string curr_string = curr_line;
        
        if (curr_string.size() > 0)
        {
            all_lines.push_back(curr_string);
        }
    }
    
    input.close();
    
    std::vector<unsigned int> lines_to_pick;
    
    unsigned int tot_size = all_lines.size();
    
    for (unsigned int i = 0; i < lines; i++)
    {
        std::cout << all_lines.at(rand()%tot_size) << "\n";
    }
    
    return 0;
}
