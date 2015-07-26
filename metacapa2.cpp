#ifdef _MSC_VER
#  pragma warning (disable : 4996) // disable -D_SCL_SECURE_NO_WARNINGS C++ 'Checked Iterators'
#endif
#include <math.h>

 #ifndef INVERT_MATRIX_HPP
 #define INVERT_MATRIX_HPP
 // REMEMBER to update "lu.hpp" header includes from boost-CVS
 #include <boost/numeric/ublas/vector.hpp>
 #include <boost/numeric/ublas/vector_proxy.hpp>
 #include <boost/numeric/ublas/matrix.hpp>
 #include <boost/numeric/ublas/triangular.hpp>
 #include <boost/numeric/ublas/lu.hpp>
 #include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas;
 /* Matrix inversion routine.
    Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
 template<class T>
bool InvertMatrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
  using namespace boost::numeric::ublas;
  typedef permutation_matrix<std::size_t> pmatrix;
  // create a working copy of the input
  matrix<T> A(input);
  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());
  // perform LU-factorization
  int res = lu_factorize(A,pm);
  if( res != 0 ) return false;
  // create identity matrix of "inverse"
  inverse.assign(ublas::identity_matrix<T>(A.size1()));
  // backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);
  return true;
}
 #endif //INVERT_MATRIX_HPP

#include <iostream>
#include <map>
#include <exception>
#include <sstream>
#include <string>
#include <fstream>

#include <boost/math/distributions/normal.hpp> // for normal_distribution.
using boost::math::normal_distribution; // default type is double.
using boost::math::normal; // typedef provides default type is double.

#include <boost/multiprecision/cpp_dec_float.hpp> // for cpp_dec_float_100
#include <boost/math/constants/constants.hpp>
#include <boost/program_options.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using std::vector;
using std::string;
using namespace std;
#include <iostream>
#include <fstream>
#include <sys/stat.h>

using namespace boost::math;
using namespace boost::multiprecision;
using namespace std;

//void info(std::string s){
//  std::cerr << "INFO:\t"<<s<<"\n";
//}

template <typename T>
void info_w(T t) 
{
  std::cerr << t << std::endl ;
}

template<typename T, typename... Args>
void info_w(T t, Args... args) // recursive variadic function
{
  std::cerr << t << " " ;

  info_w(args...) ;
}

void info_w(){
  std::cerr << endl;
}

template<typename T, typename... Args>
void info(T t, Args... args) // recursive variadic function
{
  std::cerr << "INFO:\t"<< t ;

  info_w(args...) ;
}

template <typename T>
void error_w(T t) 
{
  std::cerr << t << std::endl ;
  exit(1);
}

template<typename T, typename... Args>
void info_w(T t, Args... args) // recursive variadic function
{
  std::cerr << t << " " ;

  error_w(args...) ;
}

void error_w(){
  std::cerr << endl;
  exit(1);
}

template<typename T, typename... Args>
void error(T t, Args... args) // recursive variadic function
{
  std::cerr << "ERROR:\t"<< t ;

  error_w(args...) ;
}

long double tetrachoric(vector<bool> a, vector<bool> b){
  int n00=0;
  int n01=0;
  int n10=0;
  int n11=0;
  for(int i=0; i<a.size();i++){
    unsigned short int ap=a[i]?1:0;
    unsigned short int bp=b[i]?1:0;
    switch(2*ap-bp)
    {
      case 0:
      n00++;
      break;
      case -1:
      n01++;
      break;
      case 2:
      n10++;
      break;
      case 1:
      n11++;
    }
  }
  //DEBUG info(to_string(n00)+" "+to_string(n01)+" "+to_string(n10)+" "+to_string(n11));
  // Approximation of Digby (1973)
  long double alpha=((long double)n00*(long double)n11)/((long double)n01*(long double)n10);
  //std::cout << "n00="<<n00<<" n11="<<n11<<" n01="<<n01<<" n10="<<n10<<" alpha="<<alpha<<"\n";
  //return , boost::math::constants::pi<long double>() / 4);
  return (pow(alpha, 0.75)-1)/(pow(alpha, 0.75)+1);

}


void error(std::string s){
  std::cerr << "ERROR:\t"<<s<<"\n";
  exit(1);
}

struct position_info{
  long double beta;
  long double sdbeta;
  long double pvalue;
  bool need_precision;
  boost::multiprecision::cpp_dec_float_50 zf;
  bool bpval;
  string rsid;
  unsigned int ss;
  char allele;
};

struct pairwise{
  pairwise() : cumsum1(0),cumsum2(0), count(0),var1(0),var2(0), covar(0){}
  long double cumsum1;
  long double cumsum2;
  long int count;
  long double var1;
  long double var2;
  long double covar;
};

struct study_tuple_info {
  study_tuple_info() : initialised(false), count(0) {}
  bool initialised;
  long double var_meta_z;
  long double var_meta_beta;
  std::vector<unsigned short int> studies_map;
  unsigned short int nstudies;
  boost::numeric::ublas::vector<long double> study_weights;
  boost::numeric::ublas::vector<long double> study_weights_p;
  std::vector<std::vector<bool>> btransformed_study_vectors;
  std::vector<long double> means;
  long int count;
  boost::numeric::ublas::matrix <long double> varcovar;
  // we need a matrix to perform matrix calculations for the beta-meta analysis
  // the one above is used only to store correlations and beta variances, not for calculations.
  // the one below is an actual variance_covariance matrix, unlike the one above.
  boost::numeric::ublas::matrix <long double> varcovar_beta;

};

struct cmp_rsid{

 bool operator()(const string &a, const string &b) const
 {
  int achr=(a.at(1) == ':' ? a.at(0)- '0' : stoi(a.substr(0,2)));
  int bchr=(b.at(1) == ':' ? b.at(0)- '0' : stoi(b.substr(0,2)));
  if(a<b){return(true);}
  if(a>b){return(false);}
  long int apos;
  long int bpos;
  apos=(achr>9? stol(a.substr(3)):stol(a.substr(2)));
  bpos=(bchr>9? stol(b.substr(3)):stol(b.substr(2)));
  //  if(achr>9){info("Comparing "+a+" and "+b+" with achr="+to_string(achr)+" and bchr="+to_string(bchr)+" apos="+to_string(apos)+" bpos="+to_string(bpos));}

  return(apos < bpos);
}
};


inline std::vector<long double> produit(boost::numeric::ublas::vector<long double> &u, std::vector<long double>& v) {
  if(u.size()!=v.size()){error("Internal error 01: Mutual product of two vectors of different lengths ("+to_string(u.size())+" "+to_string(v.size())+").");}
  std::vector<long double> ret;
  for(unsigned short int i=0;i<u.size();i++){
    ret.push_back(u[i]*v[i]);
  }
  return ret;
}

inline std::vector<long double> diagonale_de(boost::numeric::ublas::matrix<long double> m){
  if(m.size1()!=m.size2()){error("Internal error 02: Attempting to extract diagonal of non-square matrix ("+to_string(m.size1())+" "+to_string(m.size2())+").");}
  std::vector<long double> ret(m.size2());
  for(unsigned short int i=0;i<m.size2();i++ ){
    ret[i]=m(i,i);
  }
  return ret;
}

inline std::vector<long double> produit(std::vector<long double> u, std::vector<long double> v) {
  if(u.size()!=v.size()){error("Internal error 01: Mutual product of two vectors of different lengths ("+to_string(u.size())+" "+to_string(v.size())+").");}
  std::vector<long double> ret;
  for(unsigned short int i=0;i<u.size();i++){
    ret.push_back(u[i]*v[i]);
  }
  return ret;
}


inline std::vector<long double> produit(boost::numeric::ublas::vector<long double> &u, boost::numeric::ublas::vector<long double> &v) {
  if(u.size()!=v.size()){error("Internal error 01: Mutual product of two vectors of different lengths ("+to_string(u.size())+" "+to_string(v.size())+").");}
  std::vector<long double> ret;
  for(unsigned short int i=0;i<u.size();i++){
    ret.push_back(u[i]*v[i]);
  }
  return ret;
}

inline std::vector<cpp_dec_float_50> produit(boost::numeric::ublas::vector<long double> &u, std::vector<cpp_dec_float_50>& v) {
  if(u.size()!=v.size()){error("Internal error 01: Mutual product of two vectors of different lengths ("+to_string(u.size())+" "+to_string(v.size())+").");}
  std::vector<cpp_dec_float_50> ret;
  for(unsigned short int i=0;i<u.size();i++){
    ret.push_back(u[i]*v[i]);
  }
  return ret;
}

inline std::vector<cpp_dec_float_50> produit(std::vector<long double> &u, std::vector<cpp_dec_float_50>& v) {
  if(u.size()!=v.size()){error("Internal error 01: Mutual product of two vectors of different lengths ("+to_string(u.size())+" "+to_string(v.size())+").");}
  std::vector<cpp_dec_float_50> ret;
  for(unsigned short int i=0;i<u.size();i++){
    ret.push_back(u[i]*v[i]);
  }
  return ret;
}

inline long double somme(std::vector<long double> u){
  long double ret=0;
  for(unsigned short int i=0;i<u.size();i++){
    ret+=u[i];
  }
  return ret;
}

inline int somme(int u[], int size){
  int ret=0;
  for(unsigned short int i=0;i<size;i++){
    ret+=u[i];
  }
  return ret;
}


inline cpp_dec_float_50 somme(std::vector<cpp_dec_float_50> u){
  cpp_dec_float_50 ret=0;
  for(unsigned short int i=0;i<u.size();i++){
    ret+=u[i];
  }
  return ret;
}


inline long double somme(boost::numeric::ublas::vector<long double> u){
  long double ret=0;
  for(unsigned short int i=0;i<u.size();i++){
    ret+=u[i];
  }
  return ret;
}

template <class T>

inline string vtostring(std::vector<T> v){
  string rval;
  for(unsigned short int i=0;i<v.size();i++){
    rval+=to_string(v[i]);
  }
  return rval;
}

inline long double produit_reciproque_asymetrique(boost::numeric::ublas::matrix<long double> m, boost::numeric::ublas::vector<long double> v){
  long double ret=0;
  for(unsigned short int i=0;i<m.size1();i++){
    for(unsigned short int j=i+1;j<m.size1(); j++){
      ret+=v[i]*v[j]*m(i,j);
    }
  }
  return ret;
}

inline long double produit_reciproque_asymetrique(boost::numeric::ublas::matrix<long double> m, std::vector<long double> v){
  long double ret=0;
  for(unsigned short int i=0;i<m.size1();i++){
    for(unsigned short int j=i+1;j<m.size1(); j++){
      ret+=v[i]*v[j]*m(i,j);
    }
  }
  return ret;
}

static void print_header(ofstream& ofs){
  ofs << "chr:pos\trsid\teffect_summary\tbeta\tbeta_se\tp-value_z_score\tz_score_sd\tp_wald_beta\tp_uncorrected\tp_province_borecki\n";
}

struct output_record
{
  string chrpos;
  long double beta;
  long double betase;
  cpp_dec_float_50 z;
  long double zse;
  cpp_dec_float_50 p;
  cpp_dec_float_50 p_wald;
  unsigned short int n;
  vector<char> status;
  vector<unsigned short int> studies;
  string rsid;
  cpp_dec_float_50 p_uncorrected;
  string info;
  void add_info(string info_a){
    if(info == ""){info=info_a;}
    else{info+=";"+info_a;}
  }
  void print(ofstream& ofs){
    // assumes that ofs is open in text append mode. If not, horrible things will happen.
    if(!ofs.is_open()){error("Internal: output file is not open.");exit(1);}
    string s(status.begin(), status.end());
    // add vtostring<unsigned short int>(studies)<< if needed
    ofs << chrpos << "\t"+rsid+"\t"+s+"\t" << to_string(beta)+"\t"+to_string(betase)+"\t"<<z<<"\t"<<zse<<"\t"<<p_wald<<"\t"<<p_uncorrected<<"\t"<<p<<"\n";
  }

};

inline bool fexists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


int main(int argc, char* argv[])
{
  using namespace std;
  namespace po = boost::program_options; 
  string FILENAME="";
  string OUTFILE="";
  vector<string> ifiles;
  po::options_description desc("METACARPA ( μ - 🐟  ): Meta-analysis in C++ Accounting for Relatedness using arbitrary Precision Arithmetic.\n===========================================================================================================\n\n\tNB: All arguments mandatory except column arguments.\n\tMATACARPA currently supports only one header line.\n\nOptions description ");
  desc.add_options()
  ("help", "This help message.")
  ("input,I", po::value<vector <string>>(), "Input file.")
  ("output,O", po::value<string>(), "Output file.")
  ("sep,t", po::value<char>(), "Input field separator. Don't forget to quote if necessary. Output field separator is always \\t.")
  ("chr-col,c", po::value<int>(), "1-based p-value column number.")
  ("pos-col,q", po::value<int>(), "1-based position column number.")
  ("all-col,a", po::value<int>(), "1-based column number for effect or reference allele.")
  ("rsid-col,r", po::value<int>(), "1-based column number for RSID or any other column that you want to keep.")
  ("pval-col,p", po::value<int>(), "1-based p-value column number.")
  ("beta-col,b", po::value<int>(), "1-based beta column number.")
  ("se-col,s", po::value<int>(), "1-based beta-SE column number.")
  ("size-col,n", po::value<int>(), "1-based sample size column (if absent, sample sizes will be assumed constant and should be appended to input file names using a comma : -I [FILENAME],[SAMPLE_SIZE]).")
  ("id-col,i", po::value<int>(), "1-based ID column number (must be unique - e.g. chr:pos-A1-A2). If absent, chr:pos will be used.")
  //    ("ss1,1", po::value<int>(), "Sample size for study 1 (in order of join).")
  //    ("ss2,2", po::value<int>(), "Sample size for study 2 (in order of join).")
  ;
  po::variables_map vm;
  //try{
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);   
  //}
  //catch (exception e){
//    error("One or more argument types were incorrect (did you forget an argument value? did you use a name where a number was needed?). Please run help with -h.");
 // }
  int sample_sizes[ifiles.size()];
  if (vm.count("help")) {
    cerr << desc << "\n";
    return 1;
  }

  if (!vm.count("input") || !vm.count("output")) {
    cerr << "ERROR:\tInput and output filename is mandatory.\n\n";
    cerr << desc << "\n";
    exit(1);
  }

  char SEP='\t';
  if(vm.count("sep")){
    SEP=vm["sep"].as<char>();
  }


  int PVAL_COL=13;
  if(vm.count("pval-col")){
    PVAL_COL=vm["pval-col"].as<int>()-1;
  }

  int BETA_COL=7;
  if(vm.count("beta-col")){
    BETA_COL=vm["beta-col"].as<int>()-1;
  }

  int SE_COL=8;
  if(vm.count("se-col")){
    SE_COL=vm["se-col"].as<int>()-1;
  }

  int SIZE_COL=-1;
  if(vm.count("size-col")){
    SIZE_COL=vm["size-col"].as<int>()-1;
  }

  int ID_COL=-1;
  if(vm.count("id-col")){
    ID_COL=vm["id-col"].as<int>()-1;
  }

  int CHR_COL=0;
  if(vm.count("chr-col")){
    CHR_COL=vm["chr-col"].as<int>()-1;
  }
 
   int POS_COL=2;
  if(vm.count("pos-col")){
    POS_COL=vm["pos-col"].as<int>()-1;
  }

  int RSID_COL=1;
  if(vm.count("id-col")){
    RSID_COL=vm["id-col"].as<int>()-1;
  }

  int ALL_COL=4;
  if(vm.count("all-col")){
    ALL_COL=vm["all-col"].as<int>()-1;
  }

  if (vm.count("input") && vm.count("output")) {
    ifiles=vm["input"].as<vector<string>>();
    if(ifiles.size()==1){error("Only one input file. Nothing to do.");}
    info("Received "+to_string(ifiles.size())+" input files:");
    for(unsigned short int i=0;i<ifiles.size();i++){
      string fn;
      stringstream filename_splitter_stream(ifiles[i]);
      getline(filename_splitter_stream, fn, ',');
      if(!fexists(fn)){error("File "+fn+" does not seem to exist.");}
      string ss;
      if(SIZE_COL==-1){
        if(filename_splitter_stream.bad()){error("METACARPA couldn't find a comma after the filename "+fn+".");}
        getline(filename_splitter_stream, ss, ',');
        int sample_size;
        try{
          sample_size=boost::lexical_cast<int>(ss);
        }catch(exception e){
          error("METACARPA doesn't think your sample size ("+ss+") is a number for file "+fn+".");
        }
        info("\t"+fn+ " with sample size "+to_string(sample_size)+".");
        sample_sizes[i]=sample_size;
      }else{info(fn," with sample size info inside.");}
      ifiles[i]=fn;

    }
    //FILENAME=vm["input"].as<string>();
    OUTFILE=vm["output"].as<string>();
    info("Writing to "+OUTFILE);
  }


  info("METACARPA ( μ - 🐟  ) successfully initialised.");
  // For each file, perform z-score calculation and binary transformation. Also update position.
  vector<int> chr_join;
  vector<long int> pos_join;
  // Initialization of the dictionary that will contain the visited positions.
  vector<map <string, position_info, cmp_rsid>> studydata;
  map <string, unsigned short int, cmp_rsid> where_present;


  // Analyse up to 15 studies
  vector<long double> means(ifiles.size());
  vector<int> meancount(ifiles.size());
  unsigned short int i=0;
  for_each(ifiles.begin(), ifiles.end(), [&](string s)-> void {
    unsigned short int current_mask=1<<i;
    string filename=s;
    ifstream inputFile(filename);
    if(!inputFile.is_open()){
      error("While opening file "+filename+", error was : "+strerror(errno)+".\n");
      exit(2);
    }
    map <string, position_info, cmp_rsid> thisstudy;
    bool firstline=true;
    info("Calculating z-scores and binomial transformation for file "+filename+"...");
    string tempLine;
    boost::math::normal_distribution<long double>  standardnormal;
    boost::math::normal_distribution<cpp_dec_float_50>  standardnormal_p;
    means[i]=0;
    int ct=0;
    while ( getline(inputFile, tempLine, '\n') ) {
      if(firstline==true){firstline=false;continue;}
      ct++;
      if(ct % 10000 == 0){cerr << ct << "\tlines\r";}
      vector <string> line;
      stringstream ss(tempLine);
      string temp;
      while (getline(ss, temp, SEP)) {  
        line.push_back(temp);
      }
      string chrpos;
      position_info thisposition;
      try{
        if(ID_COL==-1){
      chrpos=line[CHR_COL]+":"+line[POS_COL];
      }else{
        
        chrpos=line[ID_COL];
      }
    }catch(exception e){error("Something went wrong while trying to use CHR/POS or ID columns.");}
      try{thisposition.beta=stold(line[BETA_COL]);}catch(exception e){error("Something went wrong while trying to use column "+to_string(BETA_COL)+" as beta.");}
      means[i]=means[i]+thisposition.beta;
      meancount[i]++;
      try{thisposition.sdbeta=stold(line[SE_COL]);}catch(exception e){error("Something went wrong while trying to use column "+to_string(SE_COL)+" as beta-se.");}
      try{thisposition.pvalue=stold(line[PVAL_COL]);}catch(exception e){error("Something went wrong while trying to use column "+to_string(PVAL_COL)+" as p-value.");};
      try{thisposition.rsid=line[RSID_COL];}catch(exception e){error("Something went wrong while trying to use column "+to_string(RSID_COL)+" (as RSID or info).");}
//      try{thisposition.allele=boost::lexical_cast<char>(line[ALL_COL]);}catch(exception e){error("Something went wrong while trying to use column "+to_string(ALL_COL)+" as .");}
      try
      {
        thisposition.zf=(cpp_dec_float_50)(quantile(standardnormal, 1-thisposition.pvalue/2)*sign(thisposition.beta));
        thisposition.bpval=(thisposition.zf<=0);
        thisposition.need_precision=false;
      }catch(exception &ex){
        thisposition.need_precision=true;
        info("\tPosition "+chrpos+" could not be transformed with normal machine precision. Using arbitrary precision (p="+line[13]+")");
        cpp_dec_float_50 p1=1-(cpp_dec_float_50)thisposition.pvalue/2;
        thisposition.zf=quantile(standardnormal_p, p1)*sign(thisposition.beta);
        thisposition.bpval=(thisposition.zf<=0);
      }

      thisstudy[chrpos]=thisposition;
      if(where_present.count(chrpos)==0){
        where_present[chrpos]=current_mask;  
      }else{
        where_present[chrpos]=where_present[chrpos] | current_mask;
      }

    }
    studydata.push_back(thisstudy);
    means[i]=means[i]/meancount[i];
    // This is the mean beta for the whole of the file
    // It can be quite different from the mean calculated on overlapping positions only
    info("Mean effect for file "+s+" is "+to_string(means[i])+" ("+to_string(meancount[i])+" rows)");
    i++;
  });
int removed=0;


// Builds a dictionary of all available masks
map <unsigned short int, bool> mask_index;
typedef map <string, unsigned short int, cmp_rsid>::iterator it_type;
for(it_type it =where_present.begin(); it != where_present.end();){
  unsigned short int mask=it->second;
  if(!mask_index.count(mask) && (mask & (mask-1))){mask_index[mask]=true; info("Found mask "+to_string(mask));}
  ++it;
}


// Forms the study tuple data structure.
// Calculates tuple-specific binary correlations, effect means

map<unsigned short int, study_tuple_info> study_tuples_info;

info("Calculating meta-analysis specific means...");
for(auto it =where_present.begin(); it != where_present.end();++it){
  unsigned short int mask=it->second;
  string chrpos= it->first;
  for(auto it2=mask_index.begin(); it2!=mask_index.end();++it2){
    unsigned short int seen_mask=it2->first;
    if(seen_mask & mask == seen_mask){

        // If this is the first time we see this mask, we need to initialise the tuple data structure
        // Indeed it needs to know the number of studies the mask refers to before the variables can have a value
      if(!study_tuples_info[seen_mask].initialised){
        study_tuples_info[seen_mask].initialised=true;
        std::vector<unsigned short int> studies; 
        unsigned short int v=seen_mask;
        for(unsigned short int c=0;v;c++){
          if(v & 1){studies.push_back(c);}
          v>>=1;
        }
        study_tuples_info[seen_mask].studies_map=studies;
        study_tuples_info[seen_mask].nstudies=studies.size();
        study_tuples_info[seen_mask].study_weights_p=ublas::vector<long double> (studies.size(), 0);
        study_tuples_info[seen_mask].means=std::vector<long double> (study_tuples_info[seen_mask].nstudies, 0);
        study_tuples_info[seen_mask].btransformed_study_vectors=std::vector<std::vector<bool>> (study_tuples_info[seen_mask].nstudies);
        study_tuples_info[seen_mask].varcovar=boost::numeric::ublas::matrix <long double> (study_tuples_info[seen_mask].nstudies,study_tuples_info[seen_mask].nstudies,0);
        study_tuples_info[seen_mask].varcovar_beta=boost::numeric::ublas::matrix <long double> (study_tuples_info[seen_mask].nstudies,study_tuples_info[seen_mask].nstudies,0);
      }
        // update means counts
      study_tuples_info[seen_mask].count++;
      for(unsigned short int i=0; i<study_tuples_info[seen_mask].nstudies; i++){
        study_tuples_info[seen_mask].means[i]+=studydata[study_tuples_info[seen_mask].studies_map[i]][chrpos].beta;
        study_tuples_info[seen_mask].btransformed_study_vectors[i].push_back(studydata[study_tuples_info[seen_mask].studies_map[i]][chrpos].bpval);
      }

    }
  }
}

info("Calculating inter-study correlations...");
for(auto it2=mask_index.begin(); it2!=mask_index.end();++it2){
  unsigned short int seen_mask=it2->first;
  for(unsigned short int i=0; i<study_tuples_info[seen_mask].nstudies; i++){
    study_tuples_info[seen_mask].means[i]/=study_tuples_info[seen_mask].count;
    info("Mean for "+to_string(i)+"-th study within mask "+to_string(seen_mask)+" is "+to_string(study_tuples_info[seen_mask].means[i]));
    for(unsigned short int j=i+1;j<study_tuples_info[seen_mask].nstudies;j++){
      study_tuples_info[seen_mask].varcovar(i,j)=tetrachoric(study_tuples_info[seen_mask].btransformed_study_vectors[i],study_tuples_info[seen_mask].btransformed_study_vectors[j]);
      study_tuples_info[seen_mask].varcovar(j,i)=study_tuples_info[seen_mask].varcovar(i,j);
    }
  }

}

info("Calculating variance/covariance components and weights...");
for(auto it =where_present.begin(); it != where_present.end();++it){
  unsigned short int mask=it->second;
  string chrpos= it->first;
  for(auto it2=mask_index.begin(); it2!=mask_index.end();++it2){
    unsigned short int seen_mask=it2->first;
    if(seen_mask & mask == seen_mask){
      for(unsigned short int i=0; i<study_tuples_info[seen_mask].nstudies; i++){
        study_tuples_info[seen_mask].varcovar(i,i)+=pow(studydata[study_tuples_info[seen_mask].studies_map[i]][chrpos].beta-study_tuples_info[seen_mask].means[i],2);
      }
    }
  }
}



///CURRENTLY THIS STRUCTURE DOES NOT ALLOW XALCULATING WEIGHTS ON THE FLY AS NEEDED WITH A SS COLUMN

for(auto it2=mask_index.begin(); it2!=mask_index.end();++it2){
  unsigned short int seen_mask=it2->first;
  std::vector<long double> tmp_p_weights;
  for(unsigned short int i=0; i<study_tuples_info[seen_mask].nstudies; i++){
    // division step of variance calculation
    study_tuples_info[seen_mask].varcovar(i,i)/=study_tuples_info[seen_mask].count;
    // copying into beta matrix
    study_tuples_info[seen_mask].varcovar_beta(i,i)=study_tuples_info[seen_mask].varcovar(i,i);
    // weights for p-values
    //  taken from Willer et al. Bioinformatics 2010. 
    // Very surprising, because these weights don't sum to one.
    // If unsure, uncomment following line and comment Willer's.
    // study_tuples_info[seen_mask].study_weights_p.push_back(sample_sizes[i]/somme(sample_sizes, ifiles.size()));
    study_tuples_info[seen_mask].study_weights_p[i]=sqrt(sample_sizes[i])/sqrt(somme(sample_sizes, ifiles.size()));

    // this loop completes non diagonal components of the beta matrix
    for(unsigned short int j=0; j<study_tuples_info[seen_mask].nstudies;j++){
      if(i==j){continue;}
      study_tuples_info[seen_mask].varcovar_beta(i,j)=sqrt(study_tuples_info[seen_mask].varcovar_beta(i,i))*sqrt(study_tuples_info[seen_mask].varcovar_beta(j,j))*study_tuples_info[seen_mask].varcovar(i,j);
    }
  }
  // inverse of varcovar_beta matrix
  boost::numeric::ublas::matrix<long double> varcovar_beta_inverse (study_tuples_info[seen_mask].nstudies,study_tuples_info[seen_mask].nstudies);
  InvertMatrix(study_tuples_info[seen_mask].varcovar_beta,varcovar_beta_inverse);
    // weights for betas=1/sum(all terms of inverse) * line sums 
  study_tuples_info[seen_mask].study_weights=prod(boost::numeric::ublas::scalar_vector<double>(varcovar_beta_inverse.size1()), varcovar_beta_inverse);
  study_tuples_info[seen_mask].study_weights=(1/somme(study_tuples_info[seen_mask].study_weights))*study_tuples_info[seen_mask].study_weights;
    // meta beta variance
  study_tuples_info[seen_mask].var_meta_beta=somme(produit(produit(study_tuples_info[seen_mask].study_weights,study_tuples_info[seen_mask].study_weights), diagonale_de(study_tuples_info[seen_mask].varcovar_beta))) + 2*produit_reciproque_asymetrique(study_tuples_info[seen_mask].varcovar_beta, study_tuples_info[seen_mask].study_weights);
    // meta p values variance
  study_tuples_info[seen_mask].var_meta_z=somme(produit(study_tuples_info[seen_mask].study_weights_p,study_tuples_info[seen_mask].study_weights_p))+2*produit_reciproque_asymetrique(study_tuples_info[seen_mask].varcovar, study_tuples_info[seen_mask].study_weights_p);

}

info("");
info("Preliminary calculations done. Starting meta-analysis. Treating "+to_string(where_present.size())+" items.");
ofstream ofs (OUTFILE, ios::out | ios::trunc);
print_header(ofs);
// For every position ever seen, build an output record
for(auto it =where_present.begin(); it != where_present.end();++it){
  string chrpos=it->first;
  unsigned short int mask=it->second;
  output_record ord;
  ord.beta=0;
  ord.betase=0;
  ord.chrpos=chrpos;
  ord.studies=study_tuples_info[mask].studies_map;
    // Define the string of effect directions with +, - and ?
  for(unsigned short int i=0; i<ifiles.size();i++){
    if((mask & (1<<i)) == 0){
      ord.status.push_back('?');
    }else{
      //if((mask & mask-1)==0){error("mask is a power of two: "+to_string(mask)+" whereas i="+to_string(i)+" and maskcheck is "+to_string(1<<i));}
      char effect_dir=studydata[i][chrpos].beta > 0? '+' : '-';
      ord.status.push_back(effect_dir);
    }
  }

  // The output record does not need to contain anything else if the position is seen in only one study. On to the next.
  if(!(mask & (mask-1))){
    //ord.print(ofs);
    continue;}

  // meta-beta
    unsigned short int numstudies=study_tuples_info[mask].nstudies;
    std::vector<long double> beta_vector(numstudies);
    std::vector<long double> sd_vector(numstudies);
    std::vector<cpp_dec_float_50> z_vector(numstudies);
    bool need_precision=false;
    for(unsigned short k=0;k<numstudies;k++){
      beta_vector[k]=studydata[ord.studies[k]][chrpos].beta;
      z_vector[k]=studydata[ord.studies[k]][chrpos].zf;
      sd_vector[k]=studydata[ord.studies[k]][chrpos].sdbeta;
    //var_vector[k]*=var_vector[k];
    //if(isnan(var_vector[k])){
    //  info(studydata[ord.studies[k]][chrpos].sdbeta, "could not be multiplied by itself.");exit(1);}
      need_precision |= studydata[ord.studies[k]][chrpos].need_precision;
      if(ord.rsid==""){ord.rsid=studydata[ord.studies[k]][chrpos].rsid;}else{
        if(ord.rsid!=studydata[ord.studies[k]][chrpos].rsid){
          ord.rsid=ord.rsid+","+studydata[ord.studies[k]][chrpos].rsid;
        }
      }
      if(ord.rsid==""){ord.rsid="NA";}

    }
  // meta-beta
  //ord.beta=somme(produit(study_tuples_info[mask].study_weights, beta_vector));
  // meta-beta-se
  //    creating and filling the matrix.
    ublas::matrix<long double> varcovar_beta=ublas::matrix<long double>(numstudies, numstudies);
    for(unsigned short int i=0; i<numstudies; i++){
      for(unsigned short int j=0; j<numstudies;j++){
        if(i==j){varcovar_beta(i,j)=sd_vector[i]*sd_vector[i];continue;}
        varcovar_beta(i,j)=sd_vector[i]*sd_vector[j]*study_tuples_info[mask].varcovar(i,j);
      }
    }
  //    computing the weights for meta-beta analysis
  //    inverse of varcovar_beta matrix
    boost::numeric::ublas::matrix<long double> varcovar_beta_inverse (numstudies, numstudies);
    InvertMatrix(varcovar_beta,varcovar_beta_inverse);
  //    weights for betas=1/sum(all terms of inverse) * line sums 
    ublas::vector<long double> beta_weights(numstudies);
    beta_weights=prod(boost::numeric::ublas::scalar_vector<double>(varcovar_beta_inverse.size1()), varcovar_beta_inverse);
    beta_weights=(1/somme(beta_weights))*beta_weights;
  //    meta beta
    ord.beta=somme(produit(beta_weights, beta_vector));
  //    meta beta variance
    ord.betase=sqrt(somme(produit(produit(beta_weights,beta_weights), produit(sd_vector, sd_vector))) + 2*produit_reciproque_asymetrique(varcovar_beta, beta_weights));
  //
  // meta-p
  //    meta-z
    std::vector<long double> study_weights_local;    
    if(SIZE_COL==-1){
    ord.z=somme(produit(study_tuples_info[mask].study_weights_p, z_vector));
  }else{
    
    for(unsigned short int ib=0; ib<study_tuples_info[mask].nstudies; ib++){study_weights_local.push_back(sqrt(studydata[study_tuples_info[mask].studies_map[ib]][chrpos].ss));}
      //INCOMPLETE
      long double denom=sqrt(somme(produit(study_weights_local, study_weights_local)));
    for(unsigned short int ib=0; ib<study_tuples_info[mask].nstudies; ib++){study_weights_local[i]/=denom;}
      ord.z=somme(produit(study_weights_local, z_vector));


  }
  //    meta-p
    if(SIZE_COL==-1){
    ord.zse=sqrt(study_tuples_info[mask].var_meta_z);
  }else{
    ord.zse=somme(produit(study_weights_local, study_weights_local))+2*produit_reciproque_asymetrique(study_tuples_info[mask].varcovar, study_weights_local);
  }
    cpp_dec_float_50 metap;
    boost::math::normal_distribution<long double>  correctednormal(0.0, ord.zse);
    boost::math::normal_distribution<cpp_dec_float_50>  correctednormal_p(0.0, ord.zse);
    boost::math::normal_distribution<long double>  wald(0.0, 1);
    boost::math::normal_distribution<cpp_dec_float_50>  wald_p(0.0, 1);
    if(need_precision){
      ord.p=2*(1-cdf(correctednormal_p, abs(ord.z)));
      ord.p_uncorrected=2*(1-cdf(wald_p, abs(ord.z)));
      ord.p_wald=cdf(wald_p, -1*abs(ord.beta/ord.betase));
    }
    else{
      ord.p=2*(1-cdf(correctednormal, abs(static_cast<long double>(ord.z))));
      ord.p_uncorrected=2*(1-cdf(wald, abs(static_cast<long double>(ord.z))));
      ord.p_wald=2*cdf(wald, -1*abs(ord.beta/ord.betase));
      if(ord.p==0){

        info("Warning, position "+chrpos+" meta-analyses below long double precision. Analysing with multiprecision...");
        ord.p =2*(1-cdf(correctednormal_p, abs(ord.z)));
        ord.p_wald=cdf(wald_p, -1*abs(ord.beta/ord.betase));
        ord.p_uncorrected=2*(1-cdf(wald_p, abs(ord.z)));

      }
    }
    ord.print(ofs);

  }
  info("All done.");
  info("METACARPA ( μ - 🐟  ) swimming away.");
  info("Goodbye.");

}

























