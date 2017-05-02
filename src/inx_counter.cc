#include <iostream>
#include <string>

#include <cstdlib>
#include <cstring>

#include <fstream>

#include "cctable.h"

/* Possible CIGAR symbols */
const char letters[9] = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'};

/* Counter for CIGAR parsing */
unsigned int i;
/* Reference countdown */
int ref_c;
/* Read count */
unsigned int read_c;
/* Length of current CIGAR symbol */
std::string curr_num_str;
unsigned int curr_num;
/* CIGAR symbol or number */
int curr_char;
/* Deletion in read found? */
bool base_void;
/* Output the whole read? */
bool all_read;


char check;

bool dummy;
// std::vector<std::string> ref;
// std::vector<std::string> var;


/*

counter3 "$7" ""$2"" "$13" "$9" "$4" "$3" "$1

$1=NA12878   (ID campione)
$2=position
$3=VariantType (D,I,S)
$4=read name (che se non erro ti serviva solo per poi riprintarlo, dato che mi riprinterai tutta la line direi che diventa inutile)
$7=most left position
$9=cigar
$13=sequence

*/

/********************************************
 AGGIUNGERE: VIENE RICHIESTO I, D, S
 nel caso sia I o D, leggere 3 prima e 3 dopo
 ********************************************/

// CIGAR sequence unwrapped
std::string cigar_uw;


int main(int argc, char *argv[])
{
    all_read = false;

    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " [file] [ref base]\n";
        return -1;
    }

    all_read = true;

    ref_c = 0;

    CcTable<unsigned int> in_string;
    in_string.load_table(argv[1]);

    unsigned int lm_pos = (unsigned int)std::stoi(in_string.get_cell(0, 6, dummy)); // (unsigned int)atoi(argv[1]);
    unsigned int pos = (unsigned int)std::stoi(in_string.get_cell(0, 1, dummy)); // atoi(argv[2]);
    ref_c = (pos - lm_pos) + 1;

    if (ref_c <= 0)
    {
        std::ofstream output;
        std::string ref_file_name = ".";
        ref_file_name += in_string.get_cell(0, 0, dummy);
        ref_file_name += ".";
        ref_file_name += in_string.get_cell(0, 1, dummy);
        ref_file_name += ".refreads";
        std::string read_out;

        try
        {
            read_out = in_string.get_row(0, dummy);
        }
        catch (...)
        {
            std::cerr << "Error with " << argv[1] << "\n";
            return -1;
        }

        output.open(ref_file_name.c_str(), std::ofstream::out | std::ofstream::app);
        output << read_out << "\n";
        output.close();
        return 0;
    }

    bool same = false;
    // leftmost and pos are the same
    if (ref_c == 1)
    {
        same = true;
    }

    std::string read_seq = in_string.get_cell(0, 12, dummy);

    std::string check_str = in_string.get_cell(0, 2, dummy);
    check = check_str.at(0);

    std::string cs = in_string.get_cell(0, 8, dummy);

    unsigned int cig_seq_size = cs.size();
    char cig_seq[cig_seq_size];
    strcpy(cig_seq, cs.c_str());

    curr_num_str = "";
    curr_num = 0;
    curr_char = 0;

    for (i = 0; i < cig_seq_size; i++)
    {
        curr_char = (unsigned int)(cig_seq[i]);
        curr_num = 0;

        /* If curr_char is a number */
        if (curr_char >= 48 && curr_char <= 57)
        {
            curr_num_str.push_back(curr_char);
            continue;
        }
        else
        {
            curr_num = atoi(curr_num_str.c_str());

            if (same)
            {
                if (curr_char != check)
                {
                    std::ofstream output;
                    std::string ref_file_name = ".";
                    ref_file_name += in_string.get_cell(0, 0, dummy);
                    ref_file_name += ".";
                    ref_file_name += in_string.get_cell(0, 1, dummy);
                    ref_file_name += ".refreads";
                    std::string read_out;

                    try
                    {
                        read_out = in_string.get_row(0, dummy);
                    }
                    catch (...)
                    {
                        std::cerr << "Error with " << argv[1] << "\n";
                        return -1;
                    }

                    output.open(ref_file_name.c_str(), std::ofstream::out | std::ofstream::app);
                    output << read_out << "\n";
                    output.close();
                }
                else
                {
                    std::ofstream output;
                    std::string var_file_name = ".";
                    var_file_name += in_string.get_cell(0, 0, dummy);
                    var_file_name += ".varreads";
                    std::string read_out;

                    try
                    {
                        read_out = in_string.get_row(0, dummy);
                    }
                    catch (...)
                    {
                        std::cerr << "Error with " << argv[1] << "\n";
                        return -1;
                    }

                    output.open(var_file_name.c_str(), std::ofstream::out | std::ofstream::app);
                    output << read_out << "\n";
                    output.close();
                }

                return 0;
            }

            for (unsigned int i = 0; i < curr_num; i++)
            {
                cigar_uw.push_back(curr_char);
            }

            curr_num_str.clear();
            base_void = false;

            /* Match */
            if (curr_char == letters[0])
            {
                ref_c -= curr_num;
                read_c += curr_num;
            }

            /* Insertion */
            if (curr_char == letters[1])
            {
                read_c += curr_num;
            }

            /******************************************
            * num 5  6  7  8  9        10 11 12 13 14 *
            * ref A  C  C  T  T        A  T  C  A  G  *
            * rea    T  C  T  T  A  T  T     T  A  G  *
            * num    0  1  2  3  4  5  6     7  8  9  *
            * cig    X  =  =  =  I  I  X  D  X  =  =  *
            ******************************************/
            /**
                         inizio 6
                         base sul ref 13
                         base sulla read? risposta 9 -> A
                         CIG 1X3=2I1X1D1X1=

                         13 - 6 = 7

                      refc     readc

                         2 1 X 0
                         1 2 = 1
                         0 3 = 2
                        -1 4 = 3
                        -2 4 I 4
                         4 4 I 5
                         3 5 X 6
                         2 6 D 6
                         1 7 X 7
                         0 8 = 8
            */

            /* Deletion */
            if (curr_char == letters[2])
            {
                ref_c -= curr_num;
                base_void = true;
            }

            /* Skipped */
            if (curr_char == letters[3])
            {

            }

            /* Soft clip */
            if (curr_char == letters[4])
            {
                //ref_c -= curr_num;
                read_c += curr_num;
            }

            /* Hard clip */
            if (curr_char == letters[5])
            {
                // ref_c -= curr_num;
                // base_void = true;
            }

            /* Padding */
            if (curr_char == letters[6])
            {

            }

            /* Match */
            if (curr_char == letters[7])
            {
                ref_c -= curr_num;
                read_c += curr_num;
            }

            /* Mismatch */
            if (curr_char == letters[8])
            {
                ref_c -= curr_num;
                read_c += curr_num;
            }
        }

        if (ref_c < -3)
        {
            std::ofstream output;
            std::string var_file_name = ".";
            var_file_name += in_string.get_cell(0, 0, dummy);
            var_file_name += ".varreads";
            std::string ref_file_name = ".";
            ref_file_name += in_string.get_cell(0, 0, dummy);
            ref_file_name += ".";
            ref_file_name += in_string.get_cell(0, 1, dummy);
            ref_file_name += ".refreads";
            char curr_cig = 0;
            std::string read_out;

            try
            {
                read_out = in_string.get_row(0, dummy);
            }
            catch (...)
            {
                std::cerr << "Error with " << argv[1] << "\n";
                return -1;
            }

            try
            {
                curr_cig = cigar_uw.at(cigar_uw.size() + ref_c);
            }
            catch (...)
            {
                output.open(ref_file_name.c_str(), std::ofstream::out | std::ofstream::app);
                output << read_out << "\t" << pos << "\n";
                output.close();
                return 0;
            }



            if (i > 3)
            {
                try
                {
                    curr_cig = cigar_uw.at(cigar_uw.size() + ref_c);
                }
                catch (...)
                {
                    output.open(ref_file_name.c_str(), std::ofstream::out | std::ofstream::app);
                    output << read_out << "\n";
                    output.close();
                    return 0;
                }
            }



            if (base_void)
            {
                if (check == 'D')
                {
                    output.open(var_file_name.c_str(), std::ofstream::out | std::ofstream::app);
                    output << read_out << "\n";
                    output.close();
                    return 0;
                }
                output.open(ref_file_name.c_str(), std::ofstream::out | std::ofstream::app);
                output << read_out << "\n";
                output.close();
                return 0;
            }
            /*
            std::cerr << curr_cig << "\t" << check << "\n";
            if ((curr_cig == 'I' && check == 'I') || (curr_cig == 'D' && check == 'D'))
            {
                output.open(var_file_name.c_str(), std::ofstream::out | std::ofstream::app);
                output << read_out << "\n";
                output.close();
                return 0;
            }
            */
            if (check == 'I')
            {
                bool found = false;
                for (int i = 0; i < 6; i++)
                {
                    // std::cerr << i << "\\" << cigar_uw.size() << "\t" << cigar_uw.at(cigar_uw.size() - 1 - i) << "\n";
                    if (cigar_uw.at(cigar_uw.size() + ref_c + 3 - i) == 'I')
                    {
                        output.open(var_file_name.c_str(), std::ofstream::out | std::ofstream::app);
                        output << read_out << "\n";
                        output.close();
                        return 0;
                    }
                }
                output.open(ref_file_name.c_str(), std::ofstream::out | std::ofstream::app);
                output << read_out << "\n";
                output.close();
                return 0;
            }

            if (check == 'D')
            {
                bool found = false;
                for (int i = 0; i < 6; i++)
                {
                    if (cigar_uw.at(cigar_uw.size() + ref_c + 3 - i) == 'D')
                    {
                        output.open(var_file_name.c_str(), std::ofstream::out | std::ofstream::app);
                        output << read_out << "\n";
                        output.close();
                        return 0;
                    }
                }
                output.open(ref_file_name.c_str(), std::ofstream::out | std::ofstream::app);
                output << read_out << "\n";
                output.close();
                return 0;
            }

            if (curr_cig != 'I' && check == 'I')
            {
                output.open(ref_file_name.c_str(), std::ofstream::out | std::ofstream::app);
                output << read_out << "\n";
                output.close();
                return 0;
            }

            if (curr_cig != 'D' && check == 'D')
            {
                output.open(ref_file_name.c_str(), std::ofstream::out | std::ofstream::app);
                output << read_out << "\n";
                output.close();
                return 0;
            }

            if (!base_void)
            {
                if (read_seq.at(read_c - 1 + ref_c) != argv[2][0])
                {
                    output.open(var_file_name.c_str(), std::ofstream::out | std::ofstream::app);
                    output << read_out << "\n";
                    output.close();
                    return 0;
                }
                else
                {
                    output.open(ref_file_name.c_str(), std::ofstream::out | std::ofstream::app);
                    output << read_out << "\n";
                    output.close();
                    return 0;
                }

                return 0;
            }
        }
    }

    return 0;
}
