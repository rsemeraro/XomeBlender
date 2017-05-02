/******************************************************************************
    CCTable - a simple table parser
    Copyright (C) 2015  Valerio Orlandini

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/


#ifndef CCTABLE_H_
#define CCTABLE_H_


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstring>


#define LINE_MAX_LENGTH 65536


typedef struct
{
    std::string full_row;
    std::vector<std::string> cells;
    unsigned int cells_num;
} table_row_t;


// TUInt is an unsigned int type - unsigned int, unsigned short or unsigned long
template<class TUInt>
class CcTable
{
public:
    ~CcTable();

    bool load_table(const char *filename, char delimiter = '\t');

    TUInt nrow();
    TUInt ncol(TUInt row, bool &found);

    std::string get_cell(TUInt row, TUInt col, bool &found);
    std::string get_row(TUInt row, bool &found);

    bool find_first(const char *query, TUInt &row, TUInt &col, bool substring = false);
    bool find_all(const char *query, std::vector<TUInt> &rows, std::vector<TUInt> &cols, bool substring = false);

    bool set_cell(TUInt row, TUInt col, std::string content);
    bool set_row(TUInt row, std::vector<std::string> content);

private:
    std::vector<table_row_t> table_;
    char curr_line[LINE_MAX_LENGTH];
    char delimiter_;
};


template<class TUInt>
CcTable<TUInt>::~CcTable()
{
    table_.clear();
}

template<class TUInt>
bool CcTable<TUInt>::load_table(const char *filename, char delimiter)
{
    delimiter_ = delimiter;
    std::ifstream input_file;
    input_file.open(filename);
    if (!input_file.is_open())
    {
        std::cerr << "Unable to open " << filename << "\n";
        return false;
    }

    while (input_file.good())
    {
        input_file.getline(curr_line, LINE_MAX_LENGTH);

        table_row_t curr_row;

        curr_row.full_row = curr_line;
        curr_row.cells_num = 0;

        char c;
        TUInt i = 0;
        std::string curr_cell;

        while (i < curr_row.full_row.size())
        {
            c = curr_row.full_row.at(i);

            if (c == delimiter)
            {
                curr_row.cells.push_back(curr_cell);
                ++curr_row.cells_num;
                curr_cell.clear();
            }
            else
            {
                curr_cell.push_back(c);
            }

            ++i;
        }

        //if (curr_row.full_row.size() > 0)
        {
            curr_row.cells.push_back(curr_cell);
            ++curr_row.cells_num;
        }
        
        table_.push_back(curr_row);
    }
    
    input_file.close();

    return true;
}

template<class TUInt>
TUInt CcTable<TUInt>::nrow()
{
    return table_.size();
}

template<class TUInt>
TUInt CcTable<TUInt>::ncol(TUInt row, bool &found)
{
    if (row < table_.size())
    {
        found = true;
        return table_.at(row).cells_num;
    }

    found = false;
    return 0;
}

template<class TUInt>
std::string CcTable<TUInt>::get_cell(TUInt row, TUInt col, bool &found)
{
    if (row < table_.size())
    {
        if (col < table_.at(row).cells_num)
        {
            found = true;
            return table_.at(row).cells.at(col);
        }
    }

    found = false;
    return "";
}

template<class TUInt>
std::string CcTable<TUInt>::get_row(TUInt row, bool &found)
{
    if (row < table_.size())
    {
        found = true;
        return table_.at(row).full_row;
    }

    found = false;
    return "";
}

template<class TUInt>
bool CcTable<TUInt>::find_first(const char *query, TUInt &row, TUInt &col,
                                bool substring)
{
    TUInt ncol = 0;
    TUInt nrow = 0;

    for (nrow = 0; nrow < table_.size(); nrow++)
    {
        for (ncol = 0; ncol < table_.at(nrow).cells_num; ncol++)
        {
            if (substring)
            {
                if (table_.at(nrow).cells.at(ncol).find(query) != std::string::npos)
                {
                    row = nrow;
                    col = ncol;
                    return true;
                }
            }
            else
            {
                if (!strcmp(query, table_.at(nrow).cells.at(ncol).c_str()))
                {
                    row = nrow;
                    col = ncol;
                    return true;
                }
            }
        }
    }

    return false;
}

template<class TUInt>
bool CcTable<TUInt>::find_all(const char *query, std::vector<TUInt> &rows,
                              std::vector<TUInt> &cols, bool substring)
{
    TUInt ncol = 0;
    TUInt nrow = 0;
    bool found = false;

    for (nrow = 0; nrow < table_.size(); nrow++)
    {
        for (ncol = 0; ncol < table_.at(nrow).cells_num; ncol++)
        {
            if (substring)
            {
                if (table_.at(nrow).cells.at(ncol).find(query) != std::string::npos)
                {
                    rows.push_back(nrow);
                    cols.push_back(ncol);
                    return true;
                }
            }
            else if (!strcmp(query, table_.at(nrow).cells.at(ncol).c_str()))
            {
                rows.push_back(nrow);
                cols.push_back(ncol);
                found = true;
            }
        }
    }

    return found;
}

template<class TUInt>
bool CcTable<TUInt>::set_cell(TUInt row, TUInt col, std::string content)
{
    if (row < table_.size())
    {
        if (col < table_.at(row).cells_num)
        {
            table_.at(row).cells.at(col) = content;
            return true;
        }
    }

    return false;
}

template<class TUInt>
bool CcTable<TUInt>::set_row(TUInt row, std::vector<std::string> content)
{
    if (row < table_.size())
    {
        table_.at(row).cells = content;

        table_.at(row).full_row = "";

        for (TUInt i = 0; i < content.size() - 1; i++)
        {
            table_.at(row).full_row += content.at(i);
            table_.at(row).full_row += delimiter_;
        }

        table_.at(row).full_row += content.at(content.size() - 1);

        return true;
    }

    return false;
}


#endif // CCTABLE_H_
