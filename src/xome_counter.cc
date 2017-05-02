#include <iostream>
#include <string>

#include <cstdlib>
#include <cstring>


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


int main(int argc, char *argv[])
{
    all_read = false;

    if (argc < 5)
    {
        std::cerr << "Usage: " << argv[0] << " [leftmost position] [position to get] [read sequence] [CIGAR sequence] [optional: string to output if position is not '=']\n";
        return -1;
    }

    if (argc == 6)
    {
        all_read = true;
    }

    ref_c = 0;

    unsigned int lm_pos = (unsigned int)atoi(argv[1]);
    unsigned int pos = (unsigned int)atoi(argv[2]);
    ref_c = (pos - lm_pos) + 1;
    std::string read_seq = argv[3];

    std::string cs = argv[4];

    unsigned int cig_seq_size = cs.size();
    char cig_seq[cig_seq_size];
    strcpy(cig_seq, argv[4]);

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
            * num    0  1  2  3  4  5  6     7  8     *
            * cig    X  =  =  =  I  I  X  D  X  =     *
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

        if (ref_c == 0)
        {
            if (!base_void)
            {
                if (!all_read)
                {
                    std::cout << read_seq.at(read_c - 1) << "\n";
                }
                else
                {
                    if (read_seq.at(read_c - 1) != '=')
                    {
                        std::cout << argv[5] << "\n";
                    }
                }
            }
            else
            {
                if (!all_read)
                {
                    std::cout << "-\n";
                }
                else
                {
                    std::cout << argv[5] << "\n";
                }
            }
            return 0;
        }

        if (ref_c < 0)
        {
            if (!base_void)
            {
                if (!all_read)
                {
                    std::cout << read_seq.at(read_c - 1 + ref_c) << "\n";
                }
                else
                {
                    if (read_seq.at(read_c - 1 + ref_c) != '=')
                    {
                        std::cout << argv[5] << "\n";
                    }
                }
            }
            else
            {
                if (!all_read)
                {
                    std::cout << "-\n";
                }
                else
                {
                    std::cout << argv[5] << "\n";
                }
            }
            return 0;
        }
    }

    std::cerr << "An error occurred\n";

    return 0;
}
