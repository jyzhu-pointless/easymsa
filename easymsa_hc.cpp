#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<iostream>
#include<string>
#include<algorithm>
#include<vector>                          
#include<iterator> 
#include<ctime>
#include<cstdio>
#include<stack>
#include<dirent.h>
#include<fstream>
using namespace std;

void show_usage() {
	cout << "MSA by Hierarchical Clustering -- A Light-weight Multiple Sequence Alignment Tool" << endl
		<< "USAGE: msa_hc --type      [1/2] <default: 1> " << endl
		<< "                                1 - nucl; 2 - prot." << endl
		<< "              --score     [1/2/3/[path/to/customized/matrix]] <default: 1> " << endl
		<< "                          (nucl)1 - equivlant matrix; 2 - transition-transversion matrix; 3 - blast matrix." << endl
		<< "                          (prot)1 - BLOSUM62; 2 - PAM250" << endl
		<< "              --nseq      [integer, >=2] <default: 2> " << endl
		<< "              --directory [/path/to/directory/containing/sequences] " << endl
		<< "              --gap       [integer, >= 0] <default: 5>" << endl
		<< "              --output    [filename] " << endl
		<< "NOTE: sequences to be aligned should be at least 2, in FASTA format, the more the sequences, the less the sites of each." << endl;
}

//first parameter: 1-nucleic acid; 2-protein
int type = 1;//default

//second parameter: for nucleic acid: 1-equivalent matrix; 2-transition-transversion matrix; 3-blast matrix
int score = 0;//optional score matrix
char *score_cos;//customized score matrix
char name[2][21] = { { 'A','T','C','G','\0' },{'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','\0'} };
int mnucl[3][4][4] = {
	{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}},
	{{1,-5,-5,-1},{-5,1,-1,-5},{-5,-1,1,-5},{-1,-5,-5,1}},
	{{4,-5,-5,-5},{-5,4,-5,-5},{-5,-5,4,-5},{-5,-5,-5,4}}
};
//for protein: 1-BLOSUM62; 2-PAM250
int mprot[2][20][20] = {
	{{4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0},
	 {-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3},
	 {-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3},
	 {-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3},
	 {0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1},
	 {-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2},
	 {-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2},
	 {0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3},
	 {-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3},
	 {-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3},
	 {-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1},
	 {-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2},
	 {-1,-1,-2,-3,-1,0,-2,-3,-2,1,2	-1,5,0,-2,-1,-1,-1,-1,1},
	 {-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1},
	 {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2},
	 {1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2},
	 {0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0},
	 {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3},
	 {-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1},
	 {0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4},
	},
	{{2,-2,0,0,-2,0,0,1,-1,-1,-2,-1,-1,-3,1,1,1,-6,-3,0},
	 {-2,6,0,-1,-4,1,-1,-3,2,-2,-3,3,0,-4,0,0,-1,2,-4,-2},
	 {0,0,2,2,-4,1,1,0,2,-2,-3,1,-2,-3,0,1,0,-4,-2,-2},
	 {0,-1,2,4,-5,2,3,1,1,-2,-4,0,-3,-6,-1,0,0,-7,-4,-2},
	 {-2,-4,-4,-5,12,-5,-5,-3,-3,-2,-6,-5,-5,-4,-3,0,-2,-8,0,-2},
	 {0,1,1,2,-5,4,2,-1,3,-2,-2,1,-1,-5,0,-1,-1,-5,-4,-2},
	 {0,-1,1,3,-5,2,4,0,1,-2,-3,0,-2,-5,-1,0,0,-7,-4,-2},
	 {1,-3,0,1,-3,-1,0,5,-2,-3,-4,-2,-3,-5,0,1,0,-7,-5,-1},
	 {-1,2,2,1,-3,3,1,-2,6,-2,-2,0,-2,-2,0,-1,-1,-3,0,-2},
	 {-1,-2,-2,-2,-2,-2,-2,-3,-2,5,2,-2,2,1,-2,-1,0,-5,-1,4},
	 {-2,-3,-3,-4,-6,-2,-3,-4,-2,2,6,-3,4,2,-3,-3,-2,-2,-1,2},
	 {-1,3,1,0,-5,1,0,-2,0,-2,-3,5,0,-5,-1,0,0,-3,-4,-2},
	 {-1,0,-2,-3,-5,-1,-2,-3,-2,2,4,0,6,0,-2,-2,-1,-4,-2,2},
	 {-3,-4,-3,-6,-4,-5,-5,-5,-2,1,2,-5,0,9,-5,-3,-3,0,7,-1},
	 {1,0,0,-1,-3,0,-1,0,0,-2,-3,-1,-2,-5,6,1,0,-6,-5,-1},
	 {1,0,1,0,0,-1,0,1,-1,-1,-3,0,-2,-3,1,2,1,-2,-3,-1},
 	 {1,-1,0,0,-2,-1,0,0,-1,0,-2,0,-1,-3,0,1,3,-5,-3,0},
	 {-6,2,-4,-7,-8,-5,-7,-7,-3,-5,-2,-3,-4,0,-6,-2,-5,17,0,-6},
	 {-3,-4,-2,-4,0,-4,-4,-5,0,-1,-1,-4,-2,7,-5,-3,-3,0,10,-2},
	 {0,-2,-2,-2,-2,-2,-2,-1,-2,4,2,-2,2,-1,-1,-1,0,-6,-2,4},
	}
};

int matrix[20][20] = { 0 };//FINAL SCORE MATRIX
void get_scorefromfile(char* p) {
	FILE* f = fopen(p, "r");
	if (f == NULL) {
		show_usage();
		cout << "Error: Can't open the file: " << p << "!" << endl;
	}
	char r[200];//read matrix file in line
	int t[2] = { 4,20 };
	
	fgets(r, 200, f);
	int cnt = 0;
	char tmp;
	for (int i = 0;i < strlen(r);i++) { //row(resp. line) names of amino acid
		tmp = r[i];
		if (tmp >= 'A' && tmp <= 'Z') {
			name[type - 1][cnt++] = tmp;
		}
		if (cnt >= t[type - 1]) break;
	}

	int tran;//convert char into int
	int temp;
	for (int i = 0;i < t[type - 1];i++) {
		memset(r, 0, sizeof(r));
		fgets(r, 200, f);
		cnt = 0;
		tran = 0;
		for (int j = 0;j < strlen(r);j++) {
			if (r[j] == '-') { //scores can be minus
				temp = 0;j++;
				while (r[j] >='0' && r[j] <= '9') {//score may not be one-digit
					tran = (int)r[j++] - 48;
					temp = temp * 10 + tran;
				}
				matrix[i][cnt++] = -temp;
			}
			else if (r[j] >= '0' && r[j] <= '9') {
				temp = 0;
				while (r[j] >= '0' && r[j] <= '9') {
					tran = (int)r[j++] - 48;
					temp = temp * 10 + tran;
				}
				matrix[i][cnt++] = temp;
			}
			if (cnt >= t[type - 1]) break;
		}
	}
	fclose(f);
}

void get_scorematrix() {
	memset(matrix, 0, sizeof(matrix));

	// when score remains 0, users input a path
	if (score == 0) get_scorefromfile(score_cos);
	else if (type == 1) {
		for (int i = 0;i < 4;i++) {
			for (int j = 0;j < 4;j++) {
				matrix[i][j] = mnucl[score - 1][i][j];
			}
		}
	}
	else if (type == 2) {
		for (int i = 0;i < 20;i++) {
                        for (int j = 0;j < 20;j++) {
                                matrix[i][j] = mprot[score - 1][i][j];
                        }
                }
	}
/*
		char pam250[100] = "/www/wwwroot/easymsa/static/PAM250.txt";
		char blosum62[100] = "/www/wwwroot/easymsa/static/BLOSUM62.txt";
		if (score == 1) {
			get_scorefromfile(pam250);
		}
		else if (score == 2) {
			get_scorefromfile(blosum62);
*/	
}

//third parameter: the number of sequences to be aligned
int nseq = 2;//default

//fourth parameter: the dirctory containing FASTA files to be aligned
void get_seq(char *dir, string* seq) { // sequences will be in the *seq
	int len = strlen(dir);
	char* dir_copy = dir;
	if (dir[len - 1] != '/') { //dirctories put in end by '/' or not, better to unify it
		strcat(dir_copy, "/");
	}

	DIR* dirname = opendir(dir);
	if (dirname == NULL) {
		show_usage();cout << "ERROR: Can't open the dirctory: " << dir << "!" << endl;
	}
	struct dirent* ptr;
	char* filename;
	char frff[200];//first row in FASTA file
	FILE* f;
	int n = 0;

	while ((ptr = readdir(dirname)) != NULL) {
		filename = ptr->d_name;
		if (filename[0] != '.') {
			strcat(dir_copy, filename);
			f = fopen(dir_copy, "r");
			if (f == NULL) {
				show_usage();cout << "ERROR: Can't open the file: " << dir_copy << "!" << endl;
			}
			fgets(frff, 200, f);
			char cha;
			int cnt = 0;
			string tmp;
			while ((cha = fgetc(f)) != EOF) {
				if (cha >= 'A' && cha <= 'Z') {
					tmp += cha;
				}
			}
			seq[n] = tmp;
			n++;
			fclose(f);

			dir_copy[len + 1] = '\0';
			memset(frff, 0, sizeof(frff));
			memset(filename, 0, sizeof(filename));
		}
	}
	closedir(dirname);
}

//fifth parameter: penalty of gap, supposed to be positive
int gap = -5;//default

//sixth parameter: output the alignment results 
char outname[100] = { 0 };

//results of psa
class resofpsa {
public:
	string str1;//str is original
	string str2;
	string res1;//res is result after alignment
	string res2;
	int score;//score of psa
	int tag;//as a marker
};

//result of msa of every single sequence
class resofmsa {
public:
	string str;
	string res;
};

struct backtracking {
	int mark;//1-up; 2-up and left; 3-left
	int score;
};
typedef struct backtracking* path;

int grade(char a, char b);
resofpsa traceback(path** step, const int i, const int j, string str1, string str2, string res1, string res2, resofpsa resUnit);
resofpsa psa(string str1, string str2);
void getResOfPsaMatrix(string* s, int n, resofpsa** r);
int loc(string str, vector<resofmsa> ResOfMsa);
vector<resofmsa> adjust2(resofpsa tag, resofpsa tmp, vector<resofmsa> ResOfMsa);
vector<resofmsa> adjust1(resofpsa tag, vector<resofmsa> ResOfMsa, int no, int whoInRes);
bool complare(const resofpsa& a, const resofpsa& b) {
	return a.score > b.score;
}

int main(int argc, char** argv) {
	//record running time
	clock_t startTime, endTime;
	startTime = clock();

	if (argc == 1) {
		show_usage();return -1;
	}

	char* directory;
	int f_type = 0, f_score = 0, f_nseq = 0, f_output = 0, f_directory = 0, f_gap = 0;
	int optionIndex = 1;
	while (optionIndex < argc) {
		if (strncmp(argv[optionIndex], "--type", 6) == 0) {
			f_type = 1;
			char s = argv[++optionIndex][0];
			if (s == '1' || s == '2') {
				type = s - '0';
				optionIndex++;
			}
			else {
				show_usage();
				cout << "ERROR: Option 'type' must be 1 or 2." << endl;
				return -1;
			}
		}
		else if (strncmp(argv[optionIndex], "--score", 7) == 0) {
			f_score = 1;
			char s = argv[++optionIndex][0];
			if (s == '1' || s == '2' || s == '3') {
				score = s - '0';
			}
			else {
				score_cos = argv[optionIndex];
			}
			optionIndex++;
		}
		else if (strncmp(argv[optionIndex], "--nseq", 6) == 0){
			f_nseq = 1;
			nseq = atoi(argv[++optionIndex]);
			if (nseq < 2) {
				show_usage();
				cout << "ERROR: Option 'nseq' must be at least 2!" << endl;
				return -1;
			}
			optionIndex++;
		}
		else if (strncmp(argv[optionIndex], "--directory", 11) == 0) {
			f_directory = 1;
			directory = argv[++optionIndex];
			optionIndex++;
		}
		else if (strncmp(argv[optionIndex], "--gap", 5) == 0) {
			f_gap = 1;
			gap = atoi(argv[++optionIndex]);
			if (gap < 0) {
				show_usage();
				cout << "Warning: Option 'gap' is supposed to be positive!" << endl;
			}
			gap = -gap;
			optionIndex++;
		}
		else if (strncmp(argv[optionIndex], "--output", 11) == 0) {
			f_output = 1;
			int i;optionIndex++;
			for (i = 0;i < strlen(argv[optionIndex]);i++) { outname[i] = argv[optionIndex][i]; }
			outname[i] = '\0';
			optionIndex++;
		}
	}
	if (f_type == 0) cout << "Warning: Option 'type' is unset. Using default value 1." << endl;
	if (f_score == 0) {
		score = 1;cout << "Warning: Option 'score' is unset. Using default value 1." << endl;
	}
	if (f_nseq == 0) cout << "Warning: Option 'nseq' is unset. Using default value 2." << endl;
	if (f_gap == 0) cout << "Warning: Option 'gap' is unset. Using default value 5." << endl;
	if (f_directory == 0) {
		show_usage();
		return -1;
	}

	string* seq = new string[nseq];
	get_scorematrix();
	get_seq(directory, seq);

	//save psa result in ResOfPsa matrix
	resofpsa** ResOfPsa;
	ResOfPsa = new resofpsa * [nseq];
	for (int i = 0; i < nseq; i++){
		ResOfPsa[i] = new resofpsa[nseq];
		for (int j = 0; j < nseq;j++) {
			ResOfPsa[i][j] = resofpsa();
		}
	}

	//psa
	getResOfPsaMatrix(seq, nseq, ResOfPsa);

	//use different vectors to store the psa and msa results 
	int lpsa = ((nseq - 1) * nseq) / 2;
	vector<resofpsa> hierarchy_psa(0);
	vector<resofmsa> ResOfMsa(0);
	vector<resofpsa>::iterator it_psa;
	vector<resofmsa>::iterator it_msa;

	for (int i = 0;i < nseq;i++) {
		for (int j = i + 1;j < nseq;j++) {
			hierarchy_psa.push_back(ResOfPsa[i][j]);
		}
	}

	//find the highest grade
	sort(hierarchy_psa.begin(), hierarchy_psa.end(), complare);

	for (int i = 0;i < lpsa; i++) {
		if (nseq == ResOfMsa.size()) //Over when all sequences in ResOfMs
			break;

		//both str1 and str2 aren't in ResOfMsa
		if (loc(hierarchy_psa.at(i).str1, ResOfMsa) == -1 && loc(hierarchy_psa.at(i).str2, ResOfMsa) == -1) {
			if (ResOfMsa.size() == 0) { //if ResOfMsa is empty
				resofmsa s1, s2;//changing from ResOfPsa to ResOfMsa 
				s1.str = hierarchy_psa.at(i).str1;
				s1.res = hierarchy_psa.at(i).res1;
				s2.str = hierarchy_psa.at(i).str2;
				s2.res = hierarchy_psa.at(i).res2;
				ResOfMsa.push_back(s1);
				ResOfMsa.push_back(s2);
			}
			else {
				//align first sequence of ResOfMsa with hierarchy_psa.at(i).str1
				resofpsa tmp = psa(ResOfMsa.front().str, hierarchy_psa.at(i).str1);
				//combine str1, str2 with ResOfMsa
				ResOfMsa = adjust2(hierarchy_psa.at(i), tmp, ResOfMsa);
			}
		}

		//str1 in ResOfMsa£¬str2 isn't
		else if (loc(hierarchy_psa.at(i).str1, ResOfMsa) != -1 && loc(hierarchy_psa.at(i).str2, ResOfMsa) == -1) {
			int no = loc(hierarchy_psa.at(i).str1, ResOfMsa);
			ResOfMsa = adjust1(hierarchy_psa.at(i), ResOfMsa, no, 1);
		}

		//str2 in ResOfMsa£¬str1 isn't
		else if (loc(hierarchy_psa.at(i).str2, ResOfMsa) != -1 && loc(hierarchy_psa.at(i).str1, ResOfMsa) == -1) {
			int no = loc(hierarchy_psa.at(i).str2, ResOfMsa);
			ResOfMsa = adjust1(hierarchy_psa.at(i), ResOfMsa, no, 2);
		}
	}

	int seqcnter = 1;
	int le = 0;
	if (f_output == 1) {
		ofstream outfile;
		outfile.open(outname);
		outfile << "----------------------------------RESULT-----------------------------------" << endl;
		outfile << "BEFORE: " << endl;
		for (it_msa = ResOfMsa.begin();it_msa != ResOfMsa.end();it_msa++) {
			outfile << seqcnter << ".";
			le = it_msa->str.length();
			for (int i = 0;i < le;i++) {
				if (i % 70 == 0 && i != 0) {
					outfile << endl << "  ";
				}
				outfile << it_msa->str[i];
			}
			outfile << endl;
			seqcnter++;
		}
		seqcnter = 1;
		outfile << endl << "AFTER: " << endl;
		for (it_msa = ResOfMsa.begin();it_msa != ResOfMsa.end();it_msa++) {
			outfile << seqcnter << ".";
			le = it_msa->res.length();
			for (int i = 0;i < le;i++) {
				if (i % 70 == 0 && i != 0) {
					outfile << endl << "  ";
				}
				outfile << it_msa->res[i];
			}
			outfile << endl;
			seqcnter++;
		}

		endTime = clock();
		outfile << endl << "RUNNING TIME: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
		outfile << "--------------------------------------------------------------------" << endl;
		outfile.close();
	}
	else {
		cout << "----------------------------------RESULT-----------------------------------" << endl;
		cout << "BEFORE: " << endl;
		for (it_msa = ResOfMsa.begin();it_msa != ResOfMsa.end();it_msa++) {
			cout << seqcnter << ".";
			le = it_msa->str.length();
			for (int i = 0;i < le;i++) {
				if (i % 70 == 0 && i != 0) {
					cout << endl << "  ";
				}
				cout << it_msa->str[i];
			}
			cout << endl;
			seqcnter++;
		}
		seqcnter = 1;
		cout << endl << "AFTER: " << endl;
		for (it_msa = ResOfMsa.begin();it_msa != ResOfMsa.end();it_msa++) {
			cout << seqcnter << ".";
			le = it_msa->res.length();
			for (int i = 0;i < le;i++) {
				if (i % 70 == 0 && i != 0) {
					cout << endl << "  ";
				}
				cout << it_msa->res[i];
			}
			cout << endl;
			seqcnter++;
		}

		endTime = clock();
		cout << endl << "RUNNING TIME: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
		cout << "--------------------------------------------------------------------" << endl;
	}

	return 0;
}

//combine results of str1/str2-compare and str1/ResOfMsa-compare 
vector<resofmsa> adjust2(resofpsa ssc, resofpsa fsc, vector<resofmsa> ResOfMsa) {
	string S1 = ssc.res1;//aligned str1 in str1/str2-compare
	string S2 = ssc.res2;//aligned str2 in str1/str2-compare
	string F1 = fsc.res1;//aligned ResOfMsa in str1/ResOfMsa-compare
	string F2 = fsc.res2;//aligned str1 in str1/ResOfMsa-compare
	string A = ResOfMsa.front().res;
	string tmp = "";
	vector<resofmsa>::iterator it;

	int i = 0, j = 0;
	//start with combining str1 in ssc and fsc
	while (F2 != S1 && j < S1.length() && i < F2.length()) {
		if (F2[i] == S1[j]) {
			i++;j++;
		}
		else {
			if (F2[i] == '-') {
				S1.insert(j, "-");
				S2.insert(j, "-");
			}
			else if (S1[j] == '-') {
				F1.insert(i, "-");
				F2.insert(i, "-");
			}
		}
	}

	//when F2 or S1 has finished, fill the rest
	if (i == F2.length()) {
		for (int k = 0; k < S1.length() - j; k++) {
			tmp += "-";
		}
		F1 += tmp;F2 += tmp;
	}
	else if (j == S1.length()) {
		for (int k = 0; k < F2.length() - i; k++) {
			tmp += "-";
		}
		S1 += tmp;S2 += tmp;
	}
	
	//put the results of regulation into ResOfMsa
	tmp = "";
	i = 0, j = 0;
	//start with combining ResOfMsa in regulated and original ResOfMsa 
	while (A != F1 && i < A.length() && j < F1.length()) {
		if (A[i] == F1[j]) {
			i++;j++;
		}
		else {
			if (A[i] == '-') {
				F1.insert(j, "-");
				S1.insert(j, "-");
				S2.insert(j, "-");
			}
			else if (F1[j] == '-') {
				A.insert(i, "-");
				for (it = ResOfMsa.begin();it != ResOfMsa.end();it++) {
					it->res = it->res.insert(i, "-");
				}
			}
		}
	}

	if (i == A.length()) {
		for (int k = 0; k < F1.length() - j; k++) {
			tmp += "-";
		}
		A += tmp;
		for (it = ResOfMsa.begin();it != ResOfMsa.end();it++) {
			it->res = it->res + tmp;
		}
	}
	else if (j == F1.length()) {
		for (int k = 0; k < A.length() - i; k++) {
			tmp += "-";
		}
		F1 += tmp;
		S1 += tmp;
		S2 += tmp;
	}

	//After combination, str1 and str2 are put in ResOfMsa
	resofmsa s1, s2;
	s1.res = S1;
	s1.str = ssc.str1;
	s2.res = S2;
	s2.str = ssc.str2;
	ResOfMsa.push_back(s1);
	ResOfMsa.push_back(s2);
	return ResOfMsa;
}

vector<resofmsa> adjust1(resofpsa ssc, vector<resofmsa> ResOfMsa, int no, int whoInRes) {
	string A = ResOfMsa.at(no).res;
	string S1, S2;
	if (whoInRes == 1) {//str1 is in ResOfMsa
		S1 = ssc.res1;
		S2 = ssc.res2;
	}
	else { // str2 is in ResOfMsa
		S1 = ssc.res2;
		S2 = ssc.res1;
	}
	string tmp = "";
	vector<resofmsa>::iterator it;

	//the same as adjust2
	int i = 0, j = 0;
	while (A != S1 && i < A.length() && j < S1.length()) {
		if (A[i] == S1[j]) {
			i++;j++;
		}
		else {
			if (A[i] == '-') {
				S1.insert(j, "-");
				S2.insert(j, "-");
			}
			else if (S1[j] == '-') {
				A.insert(i, "-");
				for (it = ResOfMsa.begin();it != ResOfMsa.end();it++) {
					it->res = it->res.insert(i, "-");
				}
			}
		}
	}

	if (i == A.length()) { 
		for (int k = 0; k < S1.length() - j; k++) {
			tmp += "-";
		}
		A += tmp;
		for (it = ResOfMsa.begin();it != ResOfMsa.end();it++) {
			it->res = it->res + tmp;
		}
	}
	else if (j == S1.length()) {
		for (int k = 0; k < A.length() - i; k++) {
			tmp += "-";
		}
		S1 += tmp;
		S2 += tmp;
	}

	resofmsa s;
	s.res = S2;
	if (whoInRes == 1) {
		s.str = ssc.str2;
	}
	else {
		s.str = ssc.str1;
	}
	ResOfMsa.push_back(s);
	return ResOfMsa;
}

//find where str locates in ResOfMsa
int loc(string str, vector<resofmsa> ResOfMsa) {
	int i = 0;
	vector<resofmsa>::iterator it;
	//if str is in ResOfMsa, return its location
	for (it = ResOfMsa.begin();it != ResOfMsa.end();it++) {
		if (str == it->str)
			return i;
		i++;
	}
	//if str isn't in ResOfMsa, return -1
	return -1;
}

//operate psa
void getResOfPsaMatrix(string* s, int n, resofpsa** r) {
	if (n == 1){
		cout << "Error: Please input at least 2 sequences!" << endl;
	}

	for (int i = 0;i < n;i++) {
		for (int j = i + 1;j < n;j++) {
			r[i][j] = psa(s[i], s[j]);
		}
	}
}

int grade(char a, char b) {
	if (a == '-' || b == '-') {
		return gap;
	}
	else {
		int x = 0, y = 0, flag = 0;
		int t[2] = { 4,20 };//chose to use nucl or prot score matrix
		for (int k = 0;k < t[type - 1];k++) {
			if (a == name[type - 1][k]) {//find a, b in the score matrix
				x = k;flag++;if (flag >= 2)break;
			}
			if (b == name[type - 1][k]) {
				y = k;flag++;if (flag >= 2)break;
			}
		}
		return matrix[x][y];
	}
}

resofpsa traceback(path** dp, const int i, const int j, string s1, string s2, string res1, string res2, resofpsa r) {
	path tmp = dp[i][j];
	int dir = tmp->mark;//direction of tracing back
	if (r.tag != 1) { 
		if ((i == 0) && (j == 0)) { //the end of tracing back
			r.str1 = s1;r.str2 = s2;
			r.res1 = res1;r.res2 = res2;
			r.tag = 1;
			return r;
		}

		if (dir == 1) {    //gp up
			res1 = s1[i - 1] + res1;
			res2 = '-' + res2;
			r = traceback(dp, i - 1, j, s1, s2, res1, res2, r);
		}
		if (dir == 2) {    //go left and up
			res1 = s1[i - 1] + res1;
			res2 = s2[j - 1] + res2;
			r = traceback(dp, i - 1, j - 1, s1, s2, res1, res2, r);
		}
		if (dir == 3) {    //go left
			res1 = '-' + res1;
			res2 = s2[j - 1] + res2;
			r = traceback(dp, i, j - 1, s1, s2, res1, res2, r);
		}
		return r;
	}
	else {
		return r;
	}
}

resofpsa psa(string s1, string s2) {
	int m = s1.length();
	int n = s2.length();

	int v1, v2, v3, v;
	path** dp;//matrix of dynamic programming

	//initialization
	if ((dp = (path**)malloc(sizeof(path*) * (m + 1))) == NULL) {
		cout << "Error: Out of space!" << endl;
	}
	for (int i = 0; i <= m; i++) {
		if ((dp[i] = (path*)malloc(sizeof(path) * (n + 1))) == NULL) {
			cout << "Error: Out of space!" << endl;
		}
		for (int j = 0; j <= n; j++) {
			if ((dp[i][j] = (path)malloc(sizeof(backtracking))) == NULL) {
				cout << "Error: Out of space!" << endl;
			}
			dp[i][j]->mark = 0;
		}
	}

	dp[0][0]->score = 0;
	for (int i = 1; i <= m; i++) {
		dp[i][0]->score = gap * i;
		dp[i][0]->mark = 1;
	}
	for (int j = 1; j <= n; j++) {
		dp[0][j]->score = gap * j;
		dp[0][j]->mark = 3;
	}

	//N-W algorithm
	for (int i = 1; i <= m; i++) {
		for (int j = 1; j <= n; j++) {
			v1 = dp[i - 1][j]->score + gap;
			v2 = dp[i - 1][j - 1]->score + grade(s1[i - 1], s2[j - 1]);
			v3 = dp[i][j - 1]->score + gap;
			v = v1 > v2 ? v1 : v2;
			v = v > v3 ? v : v3;
			dp[i][j]->score = v;
			//find the origin
			if (v1 == v) dp[i][j]->mark = 1;
			if (v2 == v) dp[i][j]->mark = 2;
			if (v3 == v) dp[i][j]->mark = 3;
		}
	}

	//trace back
	resofpsa r;
	r.tag = 0;
	r = traceback(dp, m, n, s1, s2, "", "", r);
	r.score = dp[m][n]->score;

	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			free(dp[i][j]);
		}
		free(dp[i]);
	}
	free(dp);

	return r;
}
