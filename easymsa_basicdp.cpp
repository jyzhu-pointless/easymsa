#include<iostream>
#include<string>
#include<cstring>
#include<vector>
#include<fstream>
#include<algorithm>
#include<iomanip>
#define NUCL 4
#define PROT 20
#define PAM250 250
#define BLOSUM62 62
#define BLAST 8
#define TRTR 2
#define GLOBAL 6
#define LOCAL 5
#define NEGATIVE_INFTY -1048576
#define PREV_A 8
#define PREV_B 4
#define PREV_C 2
#define PREV_D 1
#define GAP_OPENING -5
#define GAP_EXTENDING -1
#define LINE 60

// TRTR: Transformation - transversion matrix

using namespace std;

int type, matrix, mode = 0;
int Type;
int inputFlag = 0;
int gapOpening = 0, gapExtending = 0;
int seqCnt = 0;
char seqFilepath[4][100] = {};
char seqSet[4][506] = {};
char seqAligned[4][2022] = {};
char* p_seqSet[4] = {};
char* p_seqAligned[4] = {};
int seqLen[4] = {};
int ****score;
int ****direction;
bool isPreviousAGap[4] = {0,0,0,0};

void alloc_array() {
    // allocate memory manually, or segmentation fault will occur
    int a,b,c;
    score = (int****)malloc(505 * sizeof(int***));
    direction = (int****)malloc(505 * sizeof(int***));
    for (a = 0; a < 505; a++) {
        score[a] = (int***)malloc(505 * sizeof(int**));
        direction[a] = (int***)malloc(505 * sizeof(int**));
        for (b = 0; b < 505; b++) {
            score[a][b] = (int**)malloc(505 * sizeof(int*));
            direction[a][b] = (int**)malloc(505 * sizeof(int*));
            for (c = 0; c < 505; c++) {
                score[a][b][c] = (int*)malloc(505 * sizeof(int));
                direction[a][b][c] = (int*)malloc(505 * sizeof(int));
            }
        }
    }
}

int trtr[4][4] = {{1,-5,-5,-1},{-5,1,-1,-5},{-5,-1,1,-5},{-1,-5,-5,1}};
int blast[4][4] = {{4,-5,-5,-5},{-5,4,-5,-5},{-5,-5,4,-5},{-5,-5,-5,4}};
int matFinal[20][20] = {};//FINAL SCORE MATRIX

char name[2][20] = { { 'A','T','C','G','\n' },{} };

int read_matrix(int matType) {
    ifstream matIn;
    int externalFileFlag = 0;
    if (matType == PAM250) {
        matIn.open("./static/PAM250.txt");
        externalFileFlag = 1;
    } else if (matType == BLOSUM62) {
        matIn.open("./static/BLOSUM62.txt");
        externalFileFlag = 1;
    } else if (matType == BLAST) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                matFinal[i][j] = blast[i][j];
            }
        }
    } else if (matType == TRTR) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                matFinal[i][j] = blast[i][j];
            }
        }
    }

    if (externalFileFlag && (!matIn.is_open())) {
		cout << "ERROR: Missing matrix file." << endl;
        return -1;
    }

	/*FILE* f = fopen(p, "r");
	if (f == NULL) {
		show_usage();
		cout << "Error: Can't open the file: " << p << "!" << endl;
	}*/
	char r[200];//read matrix file in line
	int t[2] = { 4,20 };
    Type = (type == NUCL) ? 1 : 2;

    matIn.getline(r, 200);
	
	int cnt = 0;
	char tmp;
	for (int i = 0;i < strlen(r);i++) { //row(resp. line) names of amino acid
		tmp = r[i];
		if (tmp >= 'A' && tmp <= 'Z') {
			name[Type - 1][cnt++] = tmp;
		}
		if (cnt >= t[Type - 1]) break;
	}

	int tran;//convert char into int
	int temp;
	for (int i = 0;i < t[Type - 1];i++) {
		memset(r, 0, sizeof(r));
		matIn.getline(r, 200);
		cnt = 0;
		tran = 0;
		for (int j = 0;j < strlen(r);j++) {
			if (r[j] == '-') { //scores can be minus
				temp = 0;j++;
				while (r[j] >='0' && r[j] <= '9') {//score may not be one-digit
					tran = (int)r[j++] - 48;
					temp = temp * 10 + tran;
				}
				matFinal[i][cnt++] = -temp;
			}
			else if (r[j] >= '0' && r[j] <= '9') {
				temp = 0;
				while (r[j] >= '0' && r[j] <= '9') {
					tran = (int)r[j++] - 48;
					temp = temp * 10 + tran;
				}
				matFinal[i][cnt++] = temp;
			}
			if (cnt >= t[Type - 1]) break;
		}
	}
	matIn.close();
    return 0;
}

int judge_gap_type(char axis, int position) {
    if (seqCnt <= 3 && axis == 'd') {
        return 0;
    }
    if (seqCnt == 2 && axis == 'c') {
        return 0;
    }
    if (isPreviousAGap[axis - 'a']) return gapExtending;
    else return gapOpening;
}

int compare_pairwise(char siteA, char siteB) {
    int siteAIndex, siteBIndex = -1;
    for (int i = 0; i < 20; i++) {
        if (name[Type - 1][i] == siteA) siteAIndex = i;
        if (name[Type - 1][i] == siteB) siteBIndex = i;
    }
    return matFinal[siteAIndex][siteBIndex];
}

int get_max_index(int* array) {
    int retIndex = 0;
    int maxCurrent = array[0];
    for (int i = 0; i < 16; i++) {
        if (array[i] > maxCurrent) {
            maxCurrent = array[i];
            retIndex = i;
        }
    }
    return retIndex;
}

int state_transition(int a, int b, int c, int d) {

    /**
     * @description: (a, b, c, d) are the position of the point in the 4-dim lattice.
     * 
     */

    int parentScore[16];
    for (int i = 0; i < 16; i++) {
        parentScore[i] = NEGATIVE_INFTY;
    }
    /*
        parentScore[16] = {   ( 0, 0, 0, 0),    ( 0, 0, 0,-1), 
                              ( 0, 0,-1, 0),    ( 0, 0,-1,-1),
                              ( 0,-1, 0, 0),    ( 0,-1, 0,-1),
                              ( 0,-1,-1, 0),    ( 0,-1,-1,-1),
                              (-1, 0, 0, 0),    (-1, 0, 0,-1), 
                              (-1, 0,-1, 0),    (-1, 0,-1,-1),
                              (-1,-1, 0, 0),    (-1,-1, 0,-1),
                              (-1,-1,-1, 0),    (-1,-1,-1,-1)    }
        Binary code: 8 (previous a), 4 (previous b), 2 (previous c), 1 (previous d)
        (0, 0, 0, 0) is useless.
     */

    parentScore[PREV_A] = (a >= 1) ? 
                          (score[a-1][b][c][d] + judge_gap_type('b',b) + judge_gap_type('c',c) + judge_gap_type('d',d)) : NEGATIVE_INFTY;
    parentScore[PREV_B] = (b >= 1) ? 
                          (score[a][b-1][c][d] + judge_gap_type('a',1) + judge_gap_type('c',c) + judge_gap_type('d',d)) : NEGATIVE_INFTY;
    parentScore[PREV_C] = (c >= 1) ? 
                          (score[a][b][c-1][d] + judge_gap_type('b',b) + judge_gap_type('a',a) + judge_gap_type('d',d)) : NEGATIVE_INFTY;
    parentScore[PREV_D] = (d >= 1) ? 
                          (score[a][b][c][d-1] + judge_gap_type('b',b) + judge_gap_type('c',c) + judge_gap_type('a',a)) : NEGATIVE_INFTY;
    parentScore[PREV_A + PREV_B] = (a >= 1 && b >= 1) ?
                                   (score[a-1][b-1][c][d] + 2*judge_gap_type('c',c) + 2*judge_gap_type('d',d) 
                                                          + compare_pairwise(seqSet[0][a],seqSet[1][b])) : NEGATIVE_INFTY;
    parentScore[PREV_A + PREV_C] = (a >= 1 && c >= 1) ?
                                   (score[a-1][b][c-1][d] + 2*judge_gap_type('b',b) + 2*judge_gap_type('d',d) 
                                                          + compare_pairwise(seqSet[0][a],seqSet[2][c])) : NEGATIVE_INFTY;
    parentScore[PREV_A + PREV_D] = (a >= 1 && d >= 1) ?
                                   (score[a-1][b][c][d-1] + 2*judge_gap_type('c',c) + 2*judge_gap_type('b',b) 
                                                          + compare_pairwise(seqSet[0][a],seqSet[3][d])) : NEGATIVE_INFTY;
    parentScore[PREV_B + PREV_C] = (b >= 1 && c >= 1) ?
                                   (score[a][b-1][c-1][d] + 2*judge_gap_type('a',a) + 2*judge_gap_type('d',d) 
                                                          + compare_pairwise(seqSet[2][c],seqSet[1][b])) : NEGATIVE_INFTY;
    parentScore[PREV_B + PREV_D] = (b >= 1 && d >= 1) ?
                                   (score[a][b-1][c][d-1] + 2*judge_gap_type('a',a) + 2*judge_gap_type('c',c) 
                                                          + compare_pairwise(seqSet[1][b],seqSet[3][d])) : NEGATIVE_INFTY;
    parentScore[PREV_C + PREV_D] = (c >= 1 && d >= 1) ?
                                   (score[a][b][c-1][d-1] + 2*judge_gap_type('a',a) + 2*judge_gap_type('b',b) 
                                                          + compare_pairwise(seqSet[2][c],seqSet[3][d])) : NEGATIVE_INFTY;
    parentScore[PREV_A + PREV_B + PREV_C] = (a >= 1 && b >= 1 && c >= 1) ?
                                            (score[a-1][b-1][c-1][d] + 3*judge_gap_type('d',d) 
                                                                     + compare_pairwise(seqSet[0][a],seqSet[1][b]) 
                                                                     + compare_pairwise(seqSet[0][a],seqSet[2][c]) 
                                                                     + compare_pairwise(seqSet[2][c],seqSet[1][b])) : NEGATIVE_INFTY;
    parentScore[PREV_A + PREV_B + PREV_D] = (a >= 1 && b >= 1 && d >= 1) ?
                                            (score[a-1][b-1][c][d-1] + 3*judge_gap_type('c',c) 
                                                                     + compare_pairwise(seqSet[0][a],seqSet[1][b]) 
                                                                     + compare_pairwise(seqSet[0][a],seqSet[3][d]) 
                                                                     + compare_pairwise(seqSet[3][d],seqSet[1][b])) : NEGATIVE_INFTY;
    parentScore[PREV_A + PREV_C + PREV_D] = (a >= 1 && c >= 1 && d >= 1) ?
                                            (score[a-1][b][c-1][d-1] + 3*judge_gap_type('b',b)
                                                                     + compare_pairwise(seqSet[0][a],seqSet[2][c]) 
                                                                     + compare_pairwise(seqSet[0][a],seqSet[3][d]) 
                                                                     + compare_pairwise(seqSet[2][c],seqSet[3][d])) : NEGATIVE_INFTY;
    parentScore[PREV_B + PREV_C + PREV_D] = (d >= 1 && b >= 1 && c >= 1) ?
                                            (score[a][b-1][c-1][d-1] + 3*judge_gap_type('a',a) 
                                                                     + compare_pairwise(seqSet[2][c],seqSet[1][b]) 
                                                                     + compare_pairwise(seqSet[3][d],seqSet[2][c]) 
                                                                     + compare_pairwise(seqSet[3][d],seqSet[1][b])) : NEGATIVE_INFTY;
    parentScore[PREV_A + PREV_B + PREV_C + PREV_D] = (a >= 1 && b >= 1 && c >= 1 && d >= 1) ?
                                                     (score[a-1][b-1][c-1][d-1] + compare_pairwise(seqSet[0][a],seqSet[1][b]) 
                                                                                + compare_pairwise(seqSet[0][a],seqSet[2][c]) 
                                                                                + compare_pairwise(seqSet[0][a],seqSet[3][d]) 
                                                                                + compare_pairwise(seqSet[1][b],seqSet[2][c]) 
                                                                                + compare_pairwise(seqSet[1][b],seqSet[3][d])
                                                                                + compare_pairwise(seqSet[2][c],seqSet[3][d])) : NEGATIVE_INFTY;

    int maxIndex = get_max_index(parentScore);
    int maxCurrentScore = parentScore[maxIndex];
    if (mode == LOCAL && maxCurrentScore < 0) {
        maxCurrentScore = 0;
    }
    if (a == 0 && b == 0 && c == 0 && d == 0) {
        maxCurrentScore = 0;
    }
    score[a][b][c][d] = maxCurrentScore;
    direction[a][b][c][d] = maxIndex;
    if (direction[a][b][c][d] % 2 == 0) {
        isPreviousAGap[3] = 1;
    }
    if ((direction[a][b][c][d] / 2) % 2 == 0) {
        isPreviousAGap[2] = 1;
    }
    if ((direction[a][b][c][d] / 4) % 2 == 0) {
        isPreviousAGap[1] = 1;
    }
    if ((direction[a][b][c][d] / 8) % 2 == 0) {
        isPreviousAGap[0] = 1;
    }


    /* 
        if (seqCnt == 2) { // 2-dim: a, b
            parentScore[PREV_A] = (a >= 1) ? (score[a][b][c][d] + judge_gap_type())
            parentScore[PREV_B] = (b >= 1) ? (score[a][b][c][d] + judge_gap_type())
        } else if (seqCnt == 3) {

        } else if (seqCnt == 4) {

        }
    */

    return 1;
}

void show_usage() {
    cout << "EasyMSA -- A Light-weight Multiple Sequence Alignment Tool (Version 0.1.0)" << endl
         << "(Basic Dynamic Programming Method)" << endl
	     << "USAGE: easymsa_basicdp --type [prot/nucl] <default: prot> " << endl
         << "                       --input [/path/to/seqA] [/path/to/seqB] [/path/to/seqC] ... " << endl
         << "                       --matrix [pam250/blosum62/blast/trtr] <default: pam250(prot), trtr(nucl)> " << endl
         << "                       --mode [global/local] <default: global>" << endl
         << "                       --gap-opening [integer, 0-100] <default: 5>" << endl
         << "                       --gap-extending [integer, 0-100] <default: 1>" << endl
         << "NOTE: Input sequences should be at least 2 and no more than 4, in FASTA format, with 3 - 500 sites. " << endl;
}

void traceback(int a, int b, int c, int d);

int main(int argc, char* argv[]) {
    
    /**
     * @description: input recognition
     */

    int optionIndex = 1;
    int seqNum = 0;
    seqSet[1][0] = '0', seqSet[2][0] = '0', seqSet[3][0] = '0', seqSet[0][0] = '0';
    

    while (optionIndex < argc) {
        if (strncmp(argv[optionIndex], "--type", 6) == 0){
            if (strcmp(argv[optionIndex + 1], "prot") == 0) {
                type = PROT;
                optionIndex += 2;
            } else if (strcmp(argv[optionIndex + 1], "nucl") == 0) {
                type = NUCL;
                optionIndex += 2;
            } else {
                show_usage();
                cout << "ERROR: Option `type' must be `prot' or `nucl'." << endl;
                return -1;
            }
        } else if (strncmp(argv[optionIndex], "--input", 7) == 0) {
            // inputFlag = 1;
            // int seqNum = 0;
            optionIndex += 1;
            while ((optionIndex < argc) && (strncmp(argv[optionIndex], "--", 2) != 0)) {   
                inputFlag = 1;

                strcpy(seqFilepath[seqNum], argv[optionIndex]);
                ifstream fin;
                ofstream fout;
                fin.open(seqFilepath[seqNum]);
                if (!fin.is_open()) {
                    show_usage();
                    cout << "ERROR: Cannot open file " << seqFilepath[seqNum] << "." << endl;
                    return -1;
                } else {
                    /*fin >> seqSet[seqNum];
                    cout << "ERROR" << endl;
                    return -1;*/
                    
                    while(fin >> (seqSet[seqNum] + 1)) {
                        if (seqSet[seqNum][0] == '>') {
                            continue; // ignore description row
                        } else if (strlen(seqSet[seqNum]) < 3) {
                            show_usage();
                            cout << "ERROR: Sequence length should be no less than 3." << endl;
                            return -1;
                        } else if (strlen(seqSet[seqNum]) > 80) {
                            show_usage();
                            cout << "ERROR: Sequence length should be no more than 500." << endl;
                            return -1;
                        } else break;
                    }
                }
                fin.close();
                if (strlen(seqSet[seqNum]) == 0) {
                    show_usage();
                    cout << "ERROR: Cannot read the sequence." << endl;
                    return -1;
                }
                optionIndex += 1;
                seqNum += 1;
                // cout << seqNum << endl;
            }
            // cout << seqNum << endl;
            if (seqNum <= 1) {
                show_usage();
                cout << "ERROR: At least 2 sequence files are required." << endl;
                return -1;
            }
            if (seqNum >= 5) {
                show_usage();
                cout << "ERROR: Sequence files must be no more than 4." << endl;
                return -1;
            }
            seqCnt = seqNum;
        } else if (strncmp(argv[optionIndex], "--matrix", 8) == 0) {
            if (type == PROT) {
                if (strcmp(argv[optionIndex + 1], "pam250") == 0) {
                    matrix = PAM250;
                    optionIndex += 2;
                } else if (strcmp(argv[optionIndex + 1], "blosum62") == 0) {
                    matrix = BLOSUM62;
                    optionIndex += 2;
                } else {
                    show_usage();
                    cout << "ERROR: Parameter " << argv[optionIndex + 1] << " on option `matrix' is illegal." << endl;
                    return -1;
                }
            } else if (type == NUCL) {
                if (strcmp(argv[optionIndex + 1], "trtr") == 0) {
                    matrix = TRTR;
                    optionIndex += 2;
                } else if (strcmp(argv[optionIndex + 1], "blast") == 0) {
                    matrix = BLAST;
                    optionIndex += 2;
                } else {
                    show_usage();
                    cout << "ERROR: Parameter " << argv[optionIndex + 1] << " on option `matrix' is illegal." << endl;
                    return -1;
                }
            }
        } else if (strncmp(argv[optionIndex], "--mode", 6) == 0) {
            if (strcmp(argv[optionIndex + 1], "global") == 0) {
                mode = GLOBAL;
                optionIndex += 2;
            } else if (strcmp(argv[optionIndex + 1], "local") == 0) {
                mode = LOCAL;
                optionIndex += 2;
            } else {
                show_usage();
                cout << "ERROR: Option `mode' must be `global' or `local'." << endl;
                return -1;
            }
        } else if (strncmp(argv[optionIndex], "--gap-opening", 14) == 0) {
            int tmpOpening;
            try {
                tmpOpening = atoi(argv[optionIndex + 1]);
            } catch (std::invalid_argument&) {
                show_usage();
                cout << "ERROR: Gap opening penalty is out of range (1-100)." << endl;
                return -1;
            } catch (std::out_of_range&) {
                show_usage();
                cout << "ERROR: Gap opening penalty is illegal." << endl;
                return -1;
            }
            gapOpening = tmpOpening;
            if (gapOpening < 1 || gapOpening > 100) {
                show_usage();
                cout << "ERROR: Gap opening penalty is out of range (1-100)." << endl;
                return -1;
            }
            optionIndex += 2;  
        } else if (strncmp(argv[optionIndex], "--gap-extending", 16) == 0) {
            int tmpExtending;
            try {
                tmpExtending = atoi(argv[optionIndex + 1]);
            } catch (std::invalid_argument&) {
                show_usage();
                cout << "ERROR: Gap extending penalty is out of range (1-100)." << endl;
                return -1;
            } catch (std::out_of_range&) {
                show_usage();
                cout << "ERROR: Gap extending penalty is illegal." << endl;
                return -1;
            }
            gapExtending = tmpExtending;
            if (gapExtending < 1 || gapExtending > 100) {
                show_usage();
                cout << "ERROR: Gap extending penalty is out of range (1-100)." << endl;
                return -1;
            }
            optionIndex += 2;
        } else {
            show_usage();
            cout << "ERROR: Option " << argv[optionIndex] << " is undefined." << endl;
        }
    }
    if (inputFlag == 0) {
        show_usage();
        cout << "ERROR: No given sequences." << endl;
        return -1;
    }
    cout << "TIP: Data loaded successfully." << endl;
    if (type == 0) {
        type = PROT;
        cout << "WARNING: Option `type' is unset. Using default value `prot'." << endl;
    }
    if (matrix == 0) {
        if (type == PROT) {
            matrix = PAM250;
            cout << "WARNING: Option `matrix' is unset. Using default value `pam250'." << endl;
        } else if (type == NUCL) {
            matrix = TRTR;
            cout << "WARNING: Option `matrix' is unset. Using default value `trtr'." << endl;
        }
    }
    if (mode == 0) {
        mode = GLOBAL;
        cout << "WARNING: Option `mode' is unset. Using default value `global'." << endl;
    }
    if (gapOpening == 0) {
        gapOpening = -5;
        cout << "WARNING: Gap opening penalty is unset. Using default value 5." << endl;
    }
    if (gapExtending == 0) {
        gapExtending = -1;
        cout << "WARNING: Gap extending penalty is unset. Using default value 1." << endl;
    }

    int read_res = 0;
    read_res = read_matrix(matrix);
    if (read_res == -1) {
        return -1;
    }

    cout << "TIP: Allocating memory, please wait ..." << endl;
    alloc_array();


    /**
     * @description: run dynamic programming
     */



    for (int a = 0; a <= strlen(seqSet[0]); a++) {
        for (int b = 0; b <= strlen(seqSet[1]); b++) {
            for (int c = 0; c <= strlen(seqSet[2]); c++) {
                for (int d = 0; d <= strlen(seqSet[3]); d++) {
                    state_transition(a,b,c,d);
                }
            }
        }
    }

    /**
     * @description: pointers initialization
     */

    for (int i = 0; i < 4; i++) {
        p_seqSet[i] = seqSet[i] + strlen(seqSet[i]) - 1;
        p_seqAligned[i] = seqAligned[i] + 2021;
    }



    /**
     * @description: traceback and print
     */

    

    traceback(strlen(seqSet[0])-1, strlen(seqSet[1])-1, strlen(seqSet[2])-1, strlen(seqSet[3])-1);
    int begin = 0;
	for (int i = 0; i < 2022; i ++) {
		if (seqAligned[0][i]) {
			begin = i;
			break;
		}
	}
    cout << "========== RESULTS ==========" << endl;
	for (int j = begin; j < 2022; j += LINE) {
		for (int i = 0; i < LINE; i ++) {
			if (j + i >= 2022) break;
			cout << seqAligned[0][j+i];
		}
		cout << endl;
		for (int i = 0; i < LINE; i ++) {
			if (j + i >= 2022) break;
			cout << seqAligned[1][j+i];
		}
		cout << endl;
        for (int i = 0; i < LINE; i ++) {
			if (j + i >= 2022) break;
			cout << seqAligned[2][j+i];
		}
		cout << endl;
        for (int i = 0; i < LINE; i ++) {
			if (j + i >= 2022) break;
			cout << seqAligned[3][j+i];
		}
		cout << endl << endl;
	}
    
	cout << "Score: " << score[strlen(seqSet[0])][strlen(seqSet[1])][strlen(seqSet[2])][strlen(seqSet[3])] << endl;
   
    return 0;
}

void traceback(int a, int b, int c, int d) {
    if (a == 0 && b == 0 && c == 0 && d == 0) {
        return;
    } else {
        for (int i = 0; i < 4; i++) {
            *p_seqAligned[i] = *p_seqSet[i];
            p_seqAligned[i]--;
            p_seqSet[i]--;
        }
        /*
            parentScore[16] = {   ( 0, 0, 0, 0),    ( 0, 0, 0,-1), 
                                ( 0, 0,-1, 0),    ( 0, 0,-1,-1),
                                ( 0,-1, 0, 0),    ( 0,-1, 0,-1),
                                ( 0,-1,-1, 0),    ( 0,-1,-1,-1),
                                (-1, 0, 0, 0),    (-1, 0, 0,-1), 
                                (-1, 0,-1, 0),    (-1, 0,-1,-1),
                                (-1,-1, 0, 0),    (-1,-1, 0,-1),
                                (-1,-1,-1, 0),    (-1,-1,-1,-1)    }
            Binary code: 8 (previous a), 4 (previous b), 2 (previous c), 1 (previous d)
            (0, 0, 0, 0) is useless.
        */
        int nextStep[4] = {1, 1, 1, 1};
        if (direction[a][b][c][d] % 2 == 0) {
            p_seqSet[3] ++;
            *(p_seqAligned[3] + 1) = '-';
            nextStep[3] = 0;
        }
        if ((direction[a][b][c][d] / 2) % 2 == 0) {
            p_seqSet[2] ++;
            *(p_seqAligned[2] + 1) = '-';
            nextStep[2] = 0;
        }
        if ((direction[a][b][c][d] / 4) % 2 == 0) {
            p_seqSet[1] ++;
            *(p_seqAligned[1] + 1) = '-';
            nextStep[1] = 0;
        }
        if ((direction[a][b][c][d] / 8) % 2 == 0) {
            p_seqSet[0] ++;
            *(p_seqAligned[0] + 1) = '-';
            nextStep[0] = 0;
        }
        //cout << "Test passed, next traceback" << a - nextStep[0] << b - nextStep[1] << c - nextStep[2] << d - nextStep[3] << endl;
        traceback(a - nextStep[0], b - nextStep[1], c - nextStep[2], d - nextStep[3]);
	}

    
}