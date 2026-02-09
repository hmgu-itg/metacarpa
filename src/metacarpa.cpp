#ifdef _MSC_VER
#  pragma warning (disable : 4996) // disable -D_SCL_SECURE_NO_WARNINGS C++ 'Checked Iterators'
#endif
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
 #ifndef INVERT_MATRIX_HPP
 #define INVERT_MATRIX_HPP
 // REMEMBER to update "lu.hpp" header includes from boost-CVS
 #include <boost/numeric/ublas/vector.hpp>
 #include <boost/numeric/ublas/vector_proxy.hpp>
 #include <boost/numeric/ublas/matrix.hpp>
 #include <boost/numeric/ublas/triangular.hpp>
 #include <boost/numeric/ublas/lu.hpp>
 #include <boost/numeric/ublas/io.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

namespace ublas = boost::numeric::ublas;
bool VERBOSE=false;
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

typedef number<cpp_dec_float<200> > cpp_dec_float_1000;
typedef number<cpp_dec_float<500> > cpp_dec_float_mid;
typedef number<cpp_dec_float<1000> > cpp_dec_float_hi;


//void info(std::string s){
//  std::cerr << "INFO:\t"<<s<<"\n";
//}

// Info and error functions
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
void error_w(T t, Args... args) // recursive variadic function
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

static bool find_alleles(std::vector<string> v, std::vector<string> u, string a, string b){
  // Finds if {a, b} is present in {{u}, {v}}
  if(v.size()!=u.size()){error("Internal error 02 (find_alleles): Vectors are of different sizes.");}
  for(unsigned short int k=0;k<v.size();k++){
    if(a==v[k] && b==u[k]){return true;}
  }
  return false;
}

inline long double somme(std::vector<long double> u){
  long double ret=0;
  for(unsigned short int i=0;i<u.size();i++){
    ret+=u[i];
  }
  return ret;
}

inline float somme(std::vector<int> u){
  float ret=0;
  for(unsigned short int i=0;i<u.size();i++){
    ret+=u[i];
  }
  return ret;
}


inline cpp_dec_float_1000 somme(std::vector<cpp_dec_float_1000> u){
  cpp_dec_float_1000 ret=0;
  for(unsigned short int i=0;i<u.size();i++){
    ret+=u[i];
  }
  return ret;
}

inline float somme(std::vector<float> u){
  float ret=0;
  for(unsigned short int i=0;i<u.size();i++){
    ret+=u[i];
  }
  return ret;
}


inline long double produit_reciproque_asymetrique(boost::numeric::ublas::matrix<long double> m, vector<long double> v){
  long double ret=0;
  for(unsigned short int i=0;i<m.size1();i++){
    for(unsigned short int j=i+1;j<m.size1(); j++){
      long double addendum=v[i]*v[j]*m(i,j);
      ret+= addendum;
    }
  }
  return ret;
}

// inline std::vector<long double> produit(std::vector<long double>& u, std::vector<long double>& v) {
//   if(u.size()!=v.size()){error("Internal error 01: Mutual product of two vectors of different lengths ("+to_string(u.size())+" "+to_string(v.size())+").");}
//   std::vector<long double> ret;
//   for(unsigned short int i=0;i<u.size();i++){
//     ret.push_back(u[i]*v[i]);
//   }
//   return ret;
// }

// template <class T, class U>
// inline std::vector<long double> produit(std::vector<T> u, std::vector<U>& v) {
//   if(u.size()!=v.size()){error("Internal error 01: Mutual product of two vectors of different lengths ("+to_string(u.size())+" "+to_string(v.size())+").");}
//   std::vector<long double> ret;
//   for(unsigned short int i=0;i<u.size();i++){

//     ret.push_back(u[i]*v[i]);
//   }
//   return ret;
// }

inline std::vector<long double> produit(std::vector<long double> u, std::vector<long double>& v) {
  if(u.size()!=v.size()){error("Internal error 01: Mutual product of two vectors of different lengths ("+to_string(u.size())+" "+to_string(v.size())+").");}
  std::vector<long double> ret;
  for(unsigned short int i=0;i<u.size();i++){
    ret.push_back(u[i]*v[i]);
  }
  return ret;
}

inline std::vector<float> produit(std::vector<float> u, std::vector<long double>& v) {
  if(u.size()!=v.size()){error("Internal error 01: Mutual product of two vectors of different lengths ("+to_string(u.size())+" "+to_string(v.size())+").");}
  std::vector<float> ret;
  for(unsigned short int i=0;i<u.size();i++){
    ret.push_back(u[i]*(float)v[i]);
  }
  return ret;
}

inline std::vector<cpp_dec_float_1000> produit(std::vector<cpp_dec_float_1000> u, std::vector<long double>& v) {
  if(u.size()!=v.size()){error("Internal error 01: Mutual product of two vectors of different lengths ("+to_string(u.size())+" "+to_string(v.size())+").");}
  std::vector<cpp_dec_float_1000> ret;
  for(unsigned short int i=0;i<u.size();i++){
    ret.push_back(u[i]*v[i]);
  }
  return ret;
}


float global_allele_frequency(vector<float> allele_frequencies, std::vector<long double> sample_sizes){
  float a=somme(produit(allele_frequencies, sample_sizes));
    float b=somme(sample_sizes);
  return(a/b);

}

inline vector<long double> colsum(boost::numeric::ublas::matrix<long double> & m){
  std::vector<long double> v;
  for (unsigned short int i=0;i<m.size1();i++){
    long double colsum=0;
    for(unsigned short int j=0; j<m.size1();j++){
      colsum+=m(j,i);
    }
    v.push_back(colsum);
    colsum=0;
  }
  return v;
}

inline bool fexists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

static bool sign_beta(long double beta){
  return beta >=0;
}

static bool b_transform(long double p, string rsid){
  boost::math::normal_distribution<long double>  standardnormal(0.0, 1.0);
  if(p==1) return true;
  try {
    return (quantile(standardnormal, 1-p) <= 0);
  }catch(exception &ex){
    info("\tPosition ",rsid," could not be transformed with normal machine precision. Using arbitrary precision (p=",p,")");
    try{
      boost::math::normal_distribution<cpp_dec_float_1000>  nd(0.0, 1.0);
      return (quantile(nd, 1-(cpp_dec_float_1000)p) <= 0);
    }catch(exception &ex2){
      try{
        boost::math::normal_distribution<cpp_dec_float_mid>  nd(0.0, 1.0);
        return (quantile(nd, 1-(cpp_dec_float_mid)p) <= 0);
      }catch(exception &ex3){
        try{
          boost::math::normal_distribution<cpp_dec_float_hi>  nd(0.0, 1.0);
          return (quantile(nd, 1-(cpp_dec_float_hi)p) <= 0);
        }catch(exception &ex4){
          // All tiers exhausted; result is unambiguous for extreme p
          return (p >= 0.5);
        }
      }
    }
  }
}

static cpp_dec_float_1000 z_transform(long double p, long double beta, string rsid){
  boost::math::normal_distribution<long double>  standardnormal(0.0, 1.0);
  int sign=beta>0?1:-1;
  if(p==1) return (cpp_dec_float_1000)1;
  try {
    return (cpp_dec_float_1000)(quantile(standardnormal, 1-p/2)*sign);
  }catch(exception &ex){
    info("\tPosition ",rsid," could not be transformed with normal machine precision. Using arbitrary precision (p=",p,")");
    try{
      boost::math::normal_distribution<cpp_dec_float_1000>  nd(0.0, 1.0);
      return (cpp_dec_float_1000)(quantile(nd, 1-(cpp_dec_float_1000)p/2)*sign);
    }catch(exception &ex2){
      try{
        boost::math::normal_distribution<cpp_dec_float_mid>  nd(0.0, 1.0);
        return (cpp_dec_float_1000)(quantile(nd, 1-(cpp_dec_float_mid)p/2)*sign);
      }catch(exception &ex3){
        try{
          boost::math::normal_distribution<cpp_dec_float_hi>  nd(0.0, 1.0);
          return (cpp_dec_float_1000)(quantile(nd, 1-(cpp_dec_float_hi)p/2)*sign);
        }catch(exception &ex4){
          info("\tPosition ",rsid," overflows all precision tiers. Capping z-score at +/-38.");
          return (cpp_dec_float_1000)(38.0)*sign;
        }
      }
    }
  }
}

static cpp_dec_float_1000 z_transform_fess(long double p, long double beta, string rsid){
  boost::math::normal_distribution<long double>  standardnormal(0.0, 1.0);
  int sign=beta>=0?1:-1;
  if(p==1) return (cpp_dec_float_1000)1;
  try {
    return (cpp_dec_float_1000)(quantile(standardnormal, p/2)*sign);
  }catch(exception &ex){
    info("\tPosition ",rsid," could not be transformed with normal machine precision. Using arbitrary precision (p=",p,")");
    try{
      boost::math::normal_distribution<cpp_dec_float_1000>  nd(0.0, 1.0);
      return (cpp_dec_float_1000)(quantile(nd, (cpp_dec_float_1000)p/2)*sign);
    }catch(exception &ex2){
      try{
        boost::math::normal_distribution<cpp_dec_float_mid>  nd(0.0, 1.0);
        return (cpp_dec_float_1000)(quantile(nd, (cpp_dec_float_mid)p/2)*sign);
      }catch(exception &ex3){
        try{
          boost::math::normal_distribution<cpp_dec_float_hi>  nd(0.0, 1.0);
          return (cpp_dec_float_1000)(quantile(nd, (cpp_dec_float_hi)p/2)*sign);
        }catch(exception &ex4){
          info("\tPosition ",rsid," overflows all precision tiers. Capping z-score at +/-38.");
          return (cpp_dec_float_1000)(-38.0)*sign;
        }
      }
    }
  }
}

static cpp_dec_float_1000 z_transform(long double p, string rsid){
  boost::math::normal_distribution<long double>  standardnormal(0.0, 1.0);
  if(p==1) return (cpp_dec_float_1000)1;
  try {
    return (cpp_dec_float_1000)(quantile(standardnormal, 1-p));
  }catch(exception &ex){
    info("\tPosition ",rsid," could not be transformed with normal machine precision. Using arbitrary precision (p=",p,")");
    try{
      boost::math::normal_distribution<cpp_dec_float_1000>  nd(0.0, 1.0);
      return (cpp_dec_float_1000)(quantile(nd, 1-(cpp_dec_float_1000)p));
    }catch(exception &ex2){
      try{
        boost::math::normal_distribution<cpp_dec_float_mid>  nd(0.0, 1.0);
        return (cpp_dec_float_1000)(quantile(nd, 1-(cpp_dec_float_mid)p));
      }catch(exception &ex3){
        try{
          boost::math::normal_distribution<cpp_dec_float_hi>  nd(0.0, 1.0);
          return (cpp_dec_float_1000)(quantile(nd, 1-(cpp_dec_float_hi)p));
        }catch(exception &ex4){
          info("\tPosition ",rsid," overflows all precision tiers.");
          return (cpp_dec_float_1000)(38.0);
        }
      }
    }
  }
}


inline std::vector<string> parse_tab(string s, char sep){
  vector <string> line;
  stringstream ss(s);
  string temp;
  while (getline(ss, temp, sep)) {  
    line.push_back(temp);
  }
  return line;
}

inline static unsigned short int bit_count(unsigned short int w){
  unsigned short int v=w;
  unsigned short int c;
  for (c = 0; v; c++)
  {
    v &= v - 1; 
  }
  return c;
}

static string parse_mask(unsigned short int mask){
  vector<string> v;
  while(mask){
    v.push_back(mask & 1? "1":"0");
    mask>>=1;
  }
  reverse(v.begin(), v.end());
  string s="[ ";
    for(string t : v){s+=(t+", ");}
      s+=" ]";
return(s);
}

template<class T>
void inline static print_vector(vector<T> v){
  string s="[\t";
    for(int c=0;c<v.size();c++){s+=boost::lexical_cast<std::string>(v[c])+"\t";}
      s+="]";
info(s);
}

// Class to store correlation matrix for betas
// Also contains serialization routines
class studies_correlation {
private:
  friend class  boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    using boost::serialization::make_nvp;
    ar & make_nvp("correlations", correlations);
  }
  
  struct count_data{
    count_data(): n00(0), n01(1), n10(1), n11(0){}
    int n00;
    int n01;
    int n10;
    int n11;
    void update_counts(bool a, bool b){
      if(VERBOSE){info("update_counts(", a," , ", b,") called");}
      unsigned short int ap=a?1:0;
      unsigned short int bp=b?1:0;
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
    long double get_alpha(){
      if(VERBOSE){info("n00=",n00,", n01=",n01,", n10=",n10,", n11=",n11);}
      long double alpha=((long double)n00*(long double)n11)/((long double)n01*(long double)n10);
      return alpha;
    }
  };

  ublas::matrix<count_data> corr_counts;
  ublas::matrix<long double> cormat;
  map <unsigned short int, ublas::matrix<long double>> correlations;
  map <unsigned short int, ublas::matrix<count_data>> count_matrix;
  string curmax;
  unsigned short int curmask;
  vector <bool> curbp;
public:
  studies_correlation(){
    curmax="";curmask=0;
  };

  void add_minimum_beta(string rsid, long double beta, unsigned short int study){
  // BEWARE this is getting complicated
  // This adds the sign of the beta, not any transformation of the p-value.
  // To use the transformed p-value use the function below.
    if(curmax==""){curmax=rsid;}
    if(curmax != rsid){
      error("INTERNAL/add_maximum called before update_minima (", curmax, " vs ",rsid,") .\n");
    }
    // update mask
    curmask|=study;
    // update vector of transformed p-values
    curbp.push_back(sign_beta(beta));
  }

  void add_minimum(string rsid, long double p, unsigned short int study){
    if(curmax==""){curmax=rsid;}
    if(curmax != rsid){
      error("INTERNAL/add_maximum called before update_minima (", curmax, " vs ",rsid,") .\n");
    }
    // update mask
    curmask|=study;
    // update vector of transformed p-values
    curbp.push_back(b_transform(p, rsid));
  }

  void update_minima(string maximum){
    curmax=maximum;
    if(curbp.size()>1){
      if(count_matrix.count(curmask)<1){add_mask(curmask);}
      update_mask(curmask, curbp);
    }
    curmask=0;
    curbp.clear();
  }

  void add_mask(unsigned short int mask){
    unsigned short int nstudies=bit_count(mask);
    cormat=ublas::matrix<long double>(nstudies,nstudies,0);
    corr_counts=ublas::matrix<count_data>(nstudies,nstudies);
    correlations[mask]=cormat;
    count_matrix[mask]=corr_counts;
  }

  void print(){
    for(auto it =correlations.begin(); it != correlations.end();++it){
      unsigned short int mask=it->first;
      ublas::matrix<long double> m=it->second;
      string line="";
      info("For mask", mask, "have matrix :");
      for(int i=0;i<m.size1();i++){
        line="|\t";
        for (int j=0; j<m.size1();j++){
          line+=to_string(m(i,j))+"\t";
        }
        line+="|";
        info(line);
      }
      
    }

  }
  void print(ublas::matrix<long double> m){
    string line;
    info("Got matrix of size ", m.size1());
    for(int i=0;i<m.size1();i++){
      line="|\t";
      for (int j=0; j<m.size1();j++){
        line+=to_string(m(i,j))+"\t";
      }
      line+="|";
      info(line);
    }

    

  }
  void write(string filename){
    std::ofstream ofs(filename);
    for(auto it =correlations.begin(); it != correlations.end();++it){
      unsigned short int mask=it->first;
      ublas::matrix<long double> m=it->second;
      ofs << mask << "\n";
      for(int i=0;i<m.size1();i++){
        for (int j=0; j<m.size1();j++){
          ofs << m(i,j)<<"\n";
        }
      }
    }
  }

  void read(string filename){
    std::ifstream ifs(filename);
    if (ifs.is_open()){
      string line;

      while(getline(ifs, line)){
        unsigned short int mask=boost::lexical_cast<unsigned short int>(line);
        add_mask(mask);
        unsigned short int nstudies=bit_count(mask);
        ublas::matrix<long double> m(nstudies,nstudies);
        for(int i=0;i<nstudies;i++){
          for (int j=0; j<nstudies;j++){
            getline(ifs, line);
            m(i,j)=stold(line);
          }
        }
        correlations[mask]=m;
      }
    }
    ifs.close();
  }


  void update_mask(unsigned short int mask, vector<bool> v){
    if(VERBOSE){info("update_mask() called with mask=",mask);}
    for(auto& msk : count_matrix){
        //info("msk.first=",msk.first, ", mask=", mask);
      if((msk.first & mask) == msk.first){
        // 1 0 0 1 1
        // 1 0 0 1 0
        //         ^
        // . . . ^
        unsigned short int encompassing=mask;
        short int encompassing_index=-1;
        short int submask_index=-1;
        vector<bool> submask_vector;

        unsigned short int submask=msk.first;
        while(encompassing && submask){
          if(encompassing&1){
            encompassing_index++;
            if(submask&1){
              submask_vector.push_back(v[encompassing_index]);

            }
          }
          encompassing>>=1;
          submask>>=1;
        }


        std::vector<unsigned short int> studies;
        unsigned short int nstudies=bit_count(msk.first);
        for(unsigned short int i=0; i<nstudies;i++){
          for(unsigned int j=i+1;j<nstudies;j++){
            try{
              count_matrix[msk.first](i,j).update_counts(submask_vector[i],submask_vector[j]);
            }catch(exception e){
              info("FAILURE with i=",i, "and j=", j,"value=",submask_vector[i] ," - ", submask_vector[j]);
              exit(1);
            }
          }
        }
      }
    }
  }


  // void update_mask(unsigned short int mask, vector<bool> v){
  //   std::vector<unsigned short int> studies;
  //   unsigned short int nstudies=bit_count(mask);
  //   cout << "mask="<<mask <<  "\n";
  //   for(unsigned short int i=0; i<nstudies;i++){
  //     for(unsigned int j=i+1;j<nstudies;j++){
  //       try{
  //         cout << "i="<<i <<" j="<< j<< " v="<< v[i] << "-"<<v[j]<<"\n";
  //         count_matrix[mask](i,j).update_counts(v[i],v[j]);
  //       }catch(exception e){
  //         info("FAILURE with i=",i, "and j=", j,"value=",v[i] ," - ", v[j]);
  //         exit(1);
  //       }
  //     }
  //   }
  // }


  void compute_correlations(){
    for(auto it =count_matrix.begin(); it != count_matrix.end();++it){
      unsigned short int mask=it->first;
      ublas::matrix<count_data> m=it->second;
      unsigned short int nstudies=bit_count(mask);
      if(VERBOSE){info("For mask ", mask, " correlations are being computed:");}
      //info(nstudies);
      //cout << "nstud=" << nstudies<<"\n";
      for(unsigned short int i=0; i<nstudies;i++){
        for(unsigned int j=i+1;j<nstudies;j++){
          long double alpha=m(i,j).get_alpha();
          if(VERBOSE){info("between study ",i,"and",j, " the alpha is ",alpha);}
          correlations[mask](i,j)=(pow(alpha, 0.75)-1)/(pow(alpha, 0.75)+1);
          if(correlations[mask](i,j)==-1){correlations[mask](i,j)=0;}
        }
      }
      //info(correlations[mask]);
    }
  }

  ublas::matrix<long double> getmat(unsigned short int mask){
    return(correlations[mask]);
  }

};

static inline string get_readable_mask(unsigned short int mask, unsigned short int numstudies) {
  string ret;
    for(unsigned short int i=0; i<numstudies;i++){
    if((mask & (1<<i)) == 0){
      ret.append(" \u274C");
    }else{
      ret.append(" \u2705");
    }
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


struct cmp_rsid{

  bool operator()(const string &a, const string &b) const
  {
    int achr=(a.at(1) == ':' ? a.at(0)- '0' : stoi(a.substr(0,2)));
    int bchr=(b.at(1) == ':' ? b.at(0)- '0' : stoi(b.substr(0,2)));
    if(achr<bchr){return(true);}
    if(achr>bchr){return(false);}
    long int apos;
    long int bpos;

    apos=(achr>9? stol(a.substr(3)):stol(a.substr(2))); 
    bpos=(bchr>9? stol(b.substr(3)):stol(b.substr(2)));
    return(apos < bpos);
  }
}rsid_comparator;

struct output_record
{
  string chrpos;
  long double beta;
  long double betase;
  cpp_dec_float_1000 z;
  cpp_dec_float_1000 z_fess;
  long double zse;
  cpp_dec_float_1000 p;
  cpp_dec_float_1000 p_wald;
  cpp_dec_float_1000 p_uncorrected;
  cpp_dec_float_1000 p_fess;
  float af;
  string a1;
  string a2;
  unsigned short int n;
  vector<char> status;
  vector<unsigned short int> studies;
  string rsid;
  int size;

  void print(ofstream& ofs){
    // assumes that ofs is open in text append mode. If not, horrible things will happen.
    if(!ofs.is_open()){error("Internal: output file is not open.");exit(1);}
    string s(status.begin(), status.end());
    ofs  << rsid <<"\t"+chrpos+"\t" <<a1<<"\t"<<a2<<"\t"<< af  << "\t"+s+"\t" << beta<<"\t"<<betase<<"\t"<<z<<"\t"<<zse<<"\t"<<p_wald<<"\t"<<p<<"\t"<<p_fess<<"\t"<<size<<"\n";
  }

};



inline void meta_analyse(std::vector<string> working_id, std::vector<string> working_rs, std::vector<string> working_a1, std::vector<string> working_a2, std::vector<long double> working_ps, std::vector<long double> working_betas, std::vector<long double> working_betase, std::vector<long double> working_weights, unsigned short int working_mask, std::vector<float> working_af,unsigned short int numstudies, studies_correlation correlations, ofstream& ofs){
  if(VERBOSE){
    cout << "Got working_id=";
    for (auto i: working_id)
      std::cout << i << ' ';
    cout << "\nGot working_rs=" ;
    for (auto i: working_rs)
      std::cout << i << ' ';
    cout << "\nGot working_a1=" ;
    for (auto i: working_a1)
      std::cout << i << ' ';
    cout << "\nGot working_a2=" ;
    for (auto i: working_a2)
      std::cout << i << ' ';
    cout << "\nGot working_ps=" ;
    for (auto i: working_ps)
      std::cout << i << ' ';
    cout << "\nGot working_betas=" ;
    for (auto i: working_betas)
      std::cout << i << ' ';
    cout << "\nGot working_betase=" ;
    for (auto i: working_betase)
      std::cout << i << ' ';
    cout << "\nGot working_weights=" ;
    for (auto i: working_weights)
      std::cout << i << ' ';
    cout << "\nGot working_AF=" ;
    for (auto i: working_af)
      std::cout << i << ' ';
    cout << "\nGot working_mask=" << working_mask << "\n";
    cout << "Got numstudies=" << numstudies << "\n";
  }

  output_record ord;
  ord.chrpos=working_rs[0];
  ord.rsid=working_id[0];
  ord.a1=working_a1[0];
  ord.a2=working_a2[0];
  if(VERBOSE){info("Calculating GAF");}
  ord.af=global_allele_frequency(working_af, working_weights);
  if(VERBOSE){info("Computing effect string");} 
  unsigned int j=0;
  for(unsigned short int i=0; i<numstudies;i++){
    if((working_mask & (1<<i)) == 0){
      ord.status.push_back('?');
    }else{
      //if((mask & mask-1)==0){error("mask is a power of two: "+to_string(mask)+" whereas i="+to_string(i)+" and maskcheck is "+to_string(1<<i));}
      char effect_dir=working_betas[j] > 0? '+' : '-';
      j++;
      ord.status.push_back(effect_dir);
    }
  }
  if(bit_count(working_mask)>1)
  {
      // compute weights
    int sum_ss=somme(working_weights);
    ord.size=sum_ss;
    int i=0;
      //for(long double l : working_weights){working_weights[i]/=sum_ss;i++;}

    for(long double l : working_weights){working_weights[i]=sqrt(working_weights[i]);working_weights[i]/=sqrt(sum_ss);i++;}

    if(VERBOSE){info("Computing meta-analysis statistic");}
      // compute transforms and sum thereof
    i=0;
    vector<cpp_dec_float_1000> zs(working_weights.size());
    vector<cpp_dec_float_1000> zs_fess(working_weights.size());
    for(long double l : working_ps){
      zs[i]=z_transform_fess(working_ps[i], working_betas[i],working_rs[i]);
      zs_fess[i]=z_transform_fess(working_ps[i], working_betas[i],working_rs[i]);
      i++;
    }
    ord.z=somme(produit(zs, working_weights));

        if(VERBOSE){info("(=",ord.z,")Computing uncorrected meta-analysis statistic.");}

    ord.z_fess=somme(produit(zs_fess, working_weights));
        if(VERBOSE){info("(=",ord.z_fess,")Computing SE.");    cout << "\nGot working_weights=" ;
    for (auto i: working_weights)
      std::cout << i << ' ';}


      // compute p-value SD
      // the first product in the line below sums to 1 by definition
      //ord.zse=sqrt(somme(produit(working_weights, working_weights))+produit_reciproque_asymetrique(correlations.getmat(working_mask), working_weights));
    ord.zse=sqrt(1+produit_reciproque_asymetrique(correlations.getmat(working_mask), working_weights));

     if(VERBOSE){ info("\n(=",ord.zse,")Computing matrix for beta weights.");}
  
      //compute matrix for beta weights
      //info("Working mask ", working_mask);
    ublas::matrix<long double> c=correlations.getmat(working_mask);
    if(c.size1()==0){error("Not enough variants with mask", get_readable_mask(working_mask, numstudies), ". Increase number of variants used to calculate the matrix, or check for allele issues (error occurred at position ",ord.chrpos,").");}
    for(i=0;i<c.size1();i++){
      for(unsigned short int j=0;j<c.size1();j++){
        if(i>j){c(i,j)=c(j,i);}else if (i==j){c(i,j)=1;}
        c(i,j)=c(i,j)*working_betase[i]*working_betase[j];
      }
    }
      // calculate weights
      //    * sum of columns/total sum of matrix
    ublas::matrix<long double> inverse(c.size1(),c.size1());

    if(!InvertMatrix(c,inverse)){correlations.print(c);error("Could not invert variance/covariance matrix.");}
         if(VERBOSE){ info("Original matrix");correlations.print(c);info("Inverted matrix");correlations.print(inverse);}

    working_weights=colsum(inverse);
    long double totsum=somme(working_weights);
    for(i=0;i<working_weights.size();i++) working_weights[i]/=totsum;

      // calculate betase
      std::vector<long double> working_var=produit(working_betase, working_betase);

    ord.betase=sqrt(somme(produit(produit(working_weights, working_weights), working_var)) + produit_reciproque_asymetrique(c, working_weights));

      // meta-beta
    ord.beta=somme(produit(working_weights, working_betas));
          // compute meta-analysis p-value
    //info("Z=",ord.z);
    //info("ZSE=", ord.zse);
    boost::math::normal_distribution<long double>  correctednormal(0.0, ord.zse);
    boost::math::normal_distribution<cpp_dec_float_1000>  correctednormal_p(0.0, ord.zse);
    boost::math::normal_distribution<long double>  wald(0.0, 1);
    boost::math::normal_distribution<cpp_dec_float_1000>  wald_p(0.0, 1);

    try{
     ord.p=2*(cdf(correctednormal, -1*abs(static_cast<long double>(ord.z))));
     ord.p_fess=2*cdf(wald, -1*abs(static_cast<long double>(ord.z_fess)));
     ord.p_wald=2*cdf(wald, -1*abs(ord.beta/ord.betase));
     if(ord.p==0){
       ord.p_wald=2*cdf(wald_p, -1*abs(ord.beta/ord.betase));
       ord.p=2*(cdf(correctednormal_p, -1*abs(ord.z)));
       ord.p_fess=2*cdf(wald_p, -1*abs(ord.z_fess));

     }
   } catch(exception e){
    ord.p=2*(cdf(correctednormal, -1*abs(ord.z.convert_to<long double>())));
    ord.p_wald=2*cdf(wald_p, -1*abs(ord.beta/ord.betase));
    ord.p_fess=2*cdf(wald_p, -1*abs(ord.z_fess));
  }

}
else {
      // SNP present in only 1 study, skip all calculations
  ord.rsid=working_id[0];
  ord.beta=working_betas[0];
  ord.betase=working_betase[0];
  ord.p=-1;
  ord.p_wald=working_ps[0];
  ord.p_fess=-1;
  ord.z=-1;
  ord.zse=-1;
  ord.p_uncorrected=-1;
  ord.a1=working_a1[0];
  ord.a2=working_a2[0];
  ord.size=working_weights[0];
}
ord.print(ofs);
}


char SEP='\t';
int PVAL_COL=13;
int BETA_COL=7;
int SE_COL=8;
int SIZE_COL=-1;
int ID_COL=-1;
int CHR_COL=0;
int POS_COL=2;
int RSID_COL=1;
int A1_COL=4;
int A2_COL=5;
int RAND_ID=0;
int AF_COL=-1;
int STEP=1;
bool SZTOPP=false;
bool USE_BETA_SIGN=false;
string MATRIX="";
string OUTFILE="";

vector<string> ifiles;
vector<int> sample_sizes;



int inline initialise(int argc, char* argv[]){
  using namespace std;
  namespace po = boost::program_options; 
  string FILENAME="";
  std::srand(std::time(0));
 RAND_ID=rand() % 1000 + 1989; // I was born that year and expect to live around 1000 years.

string prname=argv[0];

 po::options_description desc("METACARPA ( μ - 🐟  ): Meta-analysis in C++ Accounting for Relatedness using arbitrary Precision Arithmetic.\n\
===========================================================================================================\n\n\
  \tUsage : "+prname+" -I infile1[,size] [-I infile2[,size] ...] -O outfile --chr-col int --pos_col int --a1-col int\
--a2-col int --pval-col int --beta-col int --se-col int --af-col int [--size-col int] [--sep char] [--id-col int] [-m matrix_file] [-x] [-d]\n\n\
NB:\tMETACARPA currently supports only one header line in input files, which is ignored.\n\nOptions description ");
 desc.add_options()
 ("help", "This help message.")
 ("input,I", po::value<vector <string>>(), "(mandatory) Input file.")
 ("output,O", po::value<string>(), "(mandatory) Output file.")
 ("chr-col,c", po::value<int>(), "(mandatory) 1-based chromosome column number.")
 ("pos-col,q", po::value<int>(), "(mandatory) 1-based position column number.")
 ("a1-col,u", po::value<int>(), "(mandatory) 1-based column number for effect or reference allele.")
 ("a2-col,v", po::value<int>(), "(mandatory) 1-based column number for other allele.")
 //("rsid-col,r", po::value<int>(), "1-based column number for RSID or any other column that you want to keep.")
 ("pval-col,p", po::value<int>(), "(mandatory) 1-based p-value column number.")
 ("beta-col,b", po::value<int>(), "(mandatory) 1-based beta column number.")
 ("se-col,s", po::value<int>(), "(mandatory) 1-based beta-SE column number.")
 ("af-col,a", po::value<int>(), "(mandatory) 1-based effect allele frequency.")
  ("size-col,n", po::value<int>(), "1-based sample size column (if absent, sample sizes will be assumed constant and should be appended to input file names using a comma : -I infile,size).")
 ("sep,t", po::value<char>(), "Input field separator. Don't forget to quote if necessary. Output field separator is always \\t.")
 ("id-col,i", po::value<int>(), "1-based RSID column number (should be unique - e.g. chr:pos-A1-A2). If absent, chr:pos will be used.")
 //("step,e",po::value<int>(), "Consider only every n-th variant in the matrix calculation. Simulations showed that only 30k variants are necessary for accurate results.")
 ("matrix,m", po::value<string>(), "Path to a METACARPA-generated correlation matrix array.")
 ("stop,x", "Stop METACARPA after generating the matrix.")
 ("debug,d", "Toggles an extremely verbose output, for debugging purposes only.")
 ("use-beta-sign", "Use sign of beta for dichotomization in correlation matrix calculation (as described in Southam et al. 2017), instead of the default p-value transform.")


  //    ("ss1,1", po::value<int>(), "Sample size for study 1 (in order of join).")
  //    ("ss2,2", po::value<int>(), "Sample size for study 2 (in order of join).")
 ;
 po::variables_map vm;
  //try{
 po::store(po::parse_command_line(argc, argv, desc), vm);
 po::notify(vm);   
  //}
  //catch (exception e){


 if (vm.count("help")) {
  cerr << desc << "\n";
  return 1;
}

if(vm.count("debug")){
  VERBOSE=true;
}

if (vm.count("matrix")) {
  MATRIX=vm["matrix"].as<string>();
  if(!fexists(MATRIX)){error("Your matrix file "+MATRIX+" does not seem to exist.");}
  else{info("Supplied matrix file ", MATRIX);}
}

if (!vm.count("input") || !vm.count("output")) {
  cerr << "ERROR:\tInput and output filename is mandatory.\n\n";
  cerr << desc << "\n";
  exit(1);
}

if(vm.count("sep")){
  SEP=vm["sep"].as<char>();
}


if(vm.count("pval-col")){
  PVAL_COL=vm["pval-col"].as<int>()-1;
}

if(vm.count("step")){
  STEP=vm["step"].as<int>()-1;
}


if(vm.count("af-col")){
  AF_COL=vm["af-col"].as<int>()-1;
}


if(vm.count("beta-col")){
  BETA_COL=vm["beta-col"].as<int>()-1;
}

if(vm.count("se-col")){
  SE_COL=vm["se-col"].as<int>()-1;
}

if(vm.count("size-col")){
  SIZE_COL=vm["size-col"].as<int>()-1;
}

if(vm.count("id-col")){
  ID_COL=vm["id-col"].as<int>()-1;
}

if(vm.count("chr-col")){
  CHR_COL=vm["chr-col"].as<int>()-1;
}

if(vm.count("pos-col")){
  POS_COL=vm["pos-col"].as<int>()-1;
}

if(vm.count("id-col")){
  RSID_COL=vm["id-col"].as<int>()-1;
}

if(vm.count("a1-col")){
  A1_COL=vm["a1-col"].as<int>()-1;
}

if(vm.count("a2-col")){
  A2_COL=vm["a2-col"].as<int>()-1;
}


if(vm.count("stop")){
  SZTOPP=true;
}

if(vm.count("use-beta-sign")){
  USE_BETA_SIGN=true;
  info("Using sign of beta for dichotomization (Southam et al. 2017 method).");
}

if (vm.count("input") && vm.count("output")) {
  ifiles=vm["input"].as<vector<string>>();
  if(ifiles.size()==1){error("Only one input file. Nothing to do.");}
  sample_sizes=std::vector<int>(ifiles.size());
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
  info("METACARPA ( μ - 🐟  ) successfully initialised.");
}
return 0;
}

int main(int argc, char* argv[])
{
  // this performs argument parsing

  initialise(argc, argv);
  // open all files simultaneously, put them in an array
  
  std::vector<unique_ptr<ifstream>> filestreams;
  unsigned short int i=0;
  for_each(ifiles.begin(), ifiles.end(), [&](string s)-> void {
    string filename=s;
    unique_ptr<ifstream> inputFile(new ifstream(filename));

    if(!inputFile->is_open()){
      error("While opening file "+filename+", error was : "+strerror(errno)+".\n");
      exit(1);
    }
    filestreams.push_back(move(inputFile));
  });

  vector<bool> eofs(filestreams.size());
  vector<long int> counts(filestreams.size());

  std::vector<string> currentPos;
  std::vector<long double> currentPval;
  std::vector<long double> currentBeta;
  vector<long double> currentBetaSe;
  vector<string> currentRs;
  std::vector<string> currentA1;
  std::vector<string> currentA2;
  std::vector<float> currentAF;

    // Go to the first non-header line for all files and read the position and p-value
  i=0;
  for(auto& s : filestreams) {
    string tempLine;

    getline(*s, tempLine, '\n'); // gets rid of header

    getline(*s, tempLine, '\n');
    std::vector <string> line;
    try{
      line=parse_tab(tempLine, SEP);
      counts[i]=1;
      currentPos.push_back(line.at(CHR_COL)+":"+line.at(POS_COL));
      currentPval.push_back(stold(line[PVAL_COL]));
      currentBetaSe.push_back(stold(line[SE_COL]));
      currentBeta.push_back(stold(line[BETA_COL]));
      if(ID_COL>-1){currentRs.push_back(line[RSID_COL]);}else{currentRs.push_back(line.at(CHR_COL)+":"+line.at(POS_COL));}
      currentA1.push_back(line[A1_COL]);
      currentA2.push_back(line[A2_COL]);
    } catch(exception e){
      info("tempLine=\"",tempLine,"\"");
      info("chrcol=",CHR_COL," , POS_COL=",POS_COL," , PVAL_COL=",PVAL_COL," , SE_COL=",SE_COL," , BETA_COL=",BETA_COL," , RSID_COL=",RSID_COL);
      error("Impossible to parse first line of file ", ifiles[i], "\n\tCheck separator and column numbers for chromosome, position and p-value (", line.size(),")");

    }
    i++;
  }
  
  studies_correlation correlations;
  vector<string>::iterator min= min_element(currentPos.begin(),currentPos.end(),rsid_comparator);
  string minimum=*min;
  //std::copy(currentPos.begin(), currentPos.end(), std::ostream_iterator<string>(std::cout, " "));
  //  cout <<" ORG - minimum" <<minimum<<"\n";
  int poscount=0;
  
  if(MATRIX == ""){
    cout << "\n";
    info("FIRST PASS : Calculating variance-covariance matrix.");
    cout << "\n";

    while(1){
      if(VERBOSE){info("Treating variant ",poscount, ", minimum at ",minimum);}
      // Allele harmonization for first pass (needed for beta-sign method)
      string match_a1_pass1="";
      string match_a2_pass1="";
      for(unsigned short int j=0;j<currentPos.size();j++){
      // look for the minimum in the vector.
      // When found, convert p-value
      // add to data structure
      //advance file
        if(currentPos[j]==minimum){
          string tempLine;
        //info("j=",j);
        if(USE_BETA_SIGN){
          long double beta_for_sign = currentBeta[j];
          if(match_a1_pass1==""){
            match_a1_pass1=currentA1[j];
            match_a2_pass1=currentA2[j];
          } else if(currentA1[j]==match_a2_pass1 && currentA2[j]==match_a1_pass1){
            // Alleles are flipped relative to first file: negate beta
            beta_for_sign = -currentBeta[j];
            if(VERBOSE){info("First pass allele flip at ",currentPos[j]);}
          }
          correlations.add_minimum_beta(currentPos[j],beta_for_sign,pow(2,j));
        }else{
          correlations.add_minimum(currentPos[j],currentPval[j],pow(2,j));
        }
          if(!getline(*(filestreams[j]), tempLine, '\n')){currentPos[j]="30:1"; info("Read ",counts[j], " lines from file ", ifiles[j],".");eofs[j]=true;continue;}

          std::vector <string> line;
          try{
            counts[j]++;
            line=parse_tab(tempLine, SEP);
          // advance file
            currentPos[j]=line.at(CHR_COL)+":"+line.at(POS_COL);
            currentPval[j]=stold(line[PVAL_COL]);
            currentBeta[j]=stold(line[BETA_COL]);
            currentBetaSe[j]=stold(line[SE_COL]);
            if(ID_COL>-1){currentRs[j]=line[RSID_COL];}else{currentRs[j]=currentPos[j];}
            currentA1[j]=line[A1_COL];
            currentA2[j]=line[A2_COL];

          //cout << "update \t";
          //std::copy(currentPos.begin(), currentPos.end(), std::ostream_iterator<string>(std::cout, " "));
          //cout << "iterator "<<j<<"\n";

          } catch(exception e){
            error("Impossible to parse file ", ifiles[j], 
              "\n\tCheck separator and column numbers for chromosome, position and p-value ( line was ", tempLine, ", read fields : ",line.size(),")\n\t Please not that METACARPA does not currently support NA values.");
          }
        }
      }

      min= min_element(currentPos.begin(),currentPos.end(),rsid_comparator);
      if(rsid_comparator(*min, minimum)){error("At least one of the input files is unsorted: ", *min, "  <  ", minimum);}
      minimum=*min;                                                                                                                                                                                                                                                      
        if(minimum=="30:1"){break;}
      //if((poscount % STEP)==0){correlations.update_minima(minimum);}
      correlations.update_minima(minimum);
      poscount++;
      if((poscount % 10000) == 0){cout << "Processed " << poscount/1000 << "k variants\r"<<flush;}
    }
    correlations.compute_correlations();
    cout <<"\n";
    info("Writing to "+OUTFILE+"."+to_string(RAND_ID)+".matrix.txt");
    correlations.write(OUTFILE+"."+to_string(RAND_ID)+".matrix.txt");
    // info("Checking integrity of file...");
    // correlations.print();
    // studies_correlation sc;
    // sc.read(to_string(RAND_ID)+".matrix.txt");
    // sc.print();

    if(SZTOPP){
      info("Matrix has been generated. All done.");
      info("METACARPA ( μ - 🐟 ) swimming away.");
      info("Goodbye.");
      return(0);
    }

  // Now we reopen the files and conduct single-point analysis
  // First of all close and open.
  }else{
  // The else to if(MATRIX=="")
  // i.e. here we are using a precomputed matrix
    info("Reading matrix from ", MATRIX);
    try{
      correlations.read(MATRIX);
      info("Successfully read matrix ", MATRIX);

    }catch(exception e){
      error("Could not read ", MATRIX, "as a METACARPA-generated matrix.");
    }
  }




















  i=0;
  for_each(ifiles.begin(), ifiles.end(), [&](string s)-> void {
    string filename=s;
    filestreams[i]->close();
    unique_ptr<ifstream> inputFile(new ifstream(filename));
    filestreams[i].reset();
    filestreams[i]=move(inputFile);
    if(!filestreams[i]->is_open()){
      error("While opening file "+filename+", error was : "+strerror(errno)+".\n");
      exit(1);
    }
    i++;
  });

  // Skip header & read first line
  i=0;
  vector <int> weights_p (ifiles.size());
  currentPos.clear();
  currentPval.clear();
  currentRs.clear();
  currentBeta.clear();
  currentBetaSe.clear();
  currentA1.clear();
  currentA2.clear();

  string id;

  for(auto& s : filestreams) {
    string tempLine;

    getline(*s, tempLine, '\n'); // gets rid of header
    getline(*s, tempLine, '\n');
    std::vector <string> line;
    try{

      // for each input file, we store all the vital info (p-values, pos, etc) in arrays
      line=parse_tab(tempLine, SEP);
      counts[i]=1;
      currentPos.push_back(line.at(CHR_COL)+":"+line.at(POS_COL));
      currentPval.push_back(stold(line[PVAL_COL]));
      currentBetaSe.push_back(stold(line[SE_COL]));
      currentBeta.push_back(stold(line[BETA_COL]));
      if(ID_COL>-1){currentRs.push_back(line[RSID_COL]);}else{currentRs.push_back(line.at(CHR_COL)+":"+line.at(POS_COL));}
      currentA1.push_back(line[A1_COL]);
      currentA2.push_back(line[A2_COL]);
      currentAF.push_back(stof(line[AF_COL]));

      if(ID_COL>-1){id=line[RSID_COL];}else{id=line.at(CHR_COL)+":"+line.at(POS_COL);}
      // weights_p actually contain sample sizes, not weights.
      // This is because in the loop that follows, dividends are not known beforehand.
      if(SIZE_COL>-1){weights_p[i]=stoi(line[SIZE_COL]);}else{weights_p[i]=sample_sizes[i];}

    } catch(exception e){
      error("Impossible to parse first line of file ", ifiles[i], "\n\tCheck separator and column numbers for chromosome, position and p-value (", tempLine,")");
    }
    i++;
  }
  i=0;

  // Then we compute the smallest position
  min= min_element(currentPos.begin(),currentPos.end(),rsid_comparator);
  minimum=*min;
  cout << "\n";
  info("SECOND PASS : Single-point analysis.");
  cout << "\n";


  // Single-point analysis.
  ofstream ofs (OUTFILE, ios::out | ios::app);

  // Write headers
  ofs<<"rsid\tchr:pos\teffect_allele\tneffect_allele\teffect_allele_frequency\teffects\tbeta\tse\tz\tz_se\tp_wald\tp_corrected\tp_stouffer\tn\n";
  poscount=0;
  while(1){
    poscount++;
    if((poscount % 10000) == 0){cout << "Processed " << poscount/1000 << "k variants\r"<<flush;}

    vector<long double> working_ps;
    vector<long double> working_weights;
    vector<long double> working_betase;
    vector<long double> working_betas;
    std::vector<string> working_rs;
    std::vector<string> working_id;
    std::vector<string> working_a1;
    std::vector<string> working_a2;
    std::vector<float> working_af;
    unsigned int working_mask=0;

    std::vector<string> alternative_a1;
    std::vector<string> alternative_a2;

    string match_a1="";
    string match_a2="";
    bool catastrophe=false;
    for(unsigned short int j=0;j<currentPos.size();j++){
      if(currentPos[j]==minimum){
        if(match_a1==""){match_a1=currentA1[j];match_a2=currentA2[j];continue;}
        // if(alternative_a1.size()==0){alternative_a1.push_back(currentA1[j]);alternative_a2.push_back(currentA2[j]);continue;}
        if(currentA1[j] != match_a1 || currentA2[j]!=match_a2){
          if(VERBOSE){info("STOP: ", currentPos[j], minimum, match_a1, currentA1[j], match_a2, currentA2[j]);}
          if(currentA1[j]==match_a2 && currentA2[j]==match_a1){
            if(VERBOSE){info("Allele flip occurs.");}
            currentBeta[j]=-1*currentBeta[j];
            currentAF[j]=1-currentAF[j];
            currentA1[j]=match_a1;
            currentA2[j]=match_a2;
          }else{
        // there is an allele mismatch
          catastrophe=true;
         // break; <- commenting this line ensures all allele flips happen before allele mismatches are treated
        }

        // Add alternate alleles to vectors.
         
        }
      }
    }

    // ALLELE MISMATCH
    // advance all files matching the minimum
    if(catastrophe){
      // We keep a group of arrays containing all multiline records, to be meta-analysed afterwards:
      std::vector<string> dup_ids;
      std::vector<long double> dup_betas;
      std::vector<long double> dup_ses;
      std::vector<long double> dup_p;
      std::vector<unsigned short int> dup_studies;
      std::vector<int> dup_weights;
      std::vector<string> dup_a1;
      std::vector<string> dup_a2;
      std::vector<string> dup_rs;
      std::vector<float> dup_af;

      for(unsigned short int j=0;j<currentPos.size();j++){
        while(currentPos[j]==minimum){
          if(VERBOSE){info("File ",j, " ", currentRs[j], " ", currentPos[j], " ", currentA1[j], " ", currentA2[j]);}
          // Update with the current line that triggered the mismatch event:
          dup_ids.push_back(currentPos[j]+"_"+currentA1[j]+"_"+currentA2[j]);
          dup_betas.push_back(currentBeta[j]);
          dup_ses.push_back(currentBetaSe[j]);
          dup_p.push_back(currentPval[j]);
          dup_weights.push_back(weights_p[j]);
          dup_studies.push_back(j);
          dup_rs.push_back(currentRs[j]);
          dup_a1.push_back(currentA1[j]);
          dup_a2.push_back(currentA2[j]);
          dup_af.push_back(currentAF[j]);

          string tempLine;
          if(VERBOSE){ info("Reading file ",j);}
          if(!getline(*(filestreams[j]), tempLine, '\n')){currentPos[j]="30:1"; info("Read ",counts[j], "lines from file", ifiles[j],".");eofs[j]=true;continue;}
          std::vector <string> line;
          if(VERBOSE){info(tempLine);}
          try{
            counts[j]++;
            line=parse_tab(tempLine, SEP);
            currentPos[j]=line.at(CHR_COL)+":"+line.at(POS_COL);
            currentPval[j]=stold(line[PVAL_COL]);
            currentBetaSe[j]=stold(line[SE_COL]);
            currentBeta[j]=stold(line[BETA_COL]);
            if(ID_COL>-1){currentRs[j]=line[RSID_COL];}else{currentRs[j]=currentPos[j];}
            currentA1[j]=line[A1_COL];
            currentA2[j]=line[A2_COL];
            currentAF[j]=stof(line[AF_COL]);
            if(ID_COL>-1){id=line[RSID_COL];}else{id=line.at(CHR_COL)+":"+line.at(POS_COL);}
            weights_p[j]=SIZE_COL>-1?stoi(line[SIZE_COL]) : sample_sizes[j];
          } catch(exception e){
            error("Impossible to parse file ", ifiles[i], 
              "\n\tCheck separator and column numbers for chromosome, position and p-value (", line.size(),")");
          }
        }
      }
      if(VERBOSE){
        info("\n\n\n");
        for(auto lol : dup_ids){
          info(lol);
        }
        info("\n\n\n");
      }
      // At this point our arrays contain all the info on our position, independently of alleles.
      // 1. Extract separate IDs per allele
      // 2. For each of those, meta-analyse separately
    
        std::vector<string> distinct_ids(dup_ids);
        sort(distinct_ids.begin(), distinct_ids.end());
        vector<string>::iterator it=unique(distinct_ids.begin(), distinct_ids.end());
        distinct_ids.erase(it, distinct_ids.end());
        for (string &id : distinct_ids){
          if(VERBOSE){info("Distinct ID ", id);}
          // build a structure similar to the working_XXX.
          for(unsigned short int j=0;j<dup_ids.size();j++){
            if(dup_ids[j]==id){
              working_ps.push_back(dup_p[j]);
              working_weights.push_back(dup_weights[j]);
              working_rs.push_back(dup_ids[j]);
              working_betase.push_back(dup_ses[j]);
              working_betas.push_back(dup_betas[j]);
              working_id.push_back(dup_rs[j]);
              working_a1.push_back(dup_a1[j]);
              working_a2.push_back(dup_a2[j]);
              working_af.push_back(dup_af[j]);
              if(VERBOSE){info("id ", dup_ids[j]," Found in study ", dup_studies[j]);}
              working_mask|=(unsigned short int)(pow(2,dup_studies[j]));
            }
          }
          //meta-analyse here:
          meta_analyse(working_id, working_rs, working_a1, working_a2, working_ps, working_betas,  working_betase, working_weights, working_mask, working_af, ifiles.size(), correlations, ofs);



          working_a2.clear();
          working_a1.clear();
          working_weights.clear();
          working_id.clear();
          working_betas.clear();
          working_betase.clear();
          working_ps.clear();
          working_rs.clear();
          working_af.clear();
          working_mask=0;
        }
        min= min_element(currentPos.begin(),currentPos.end(),rsid_comparator);
        if(rsid_comparator(*min, minimum)){error("At least one of the input files is unsorted: ", *min, "  <  ", minimum);}
        minimum=*min;                                                                                                                                                                                                                                                      
        if(minimum=="30:1"){break;}
        correlations.update_minima(minimum);
        continue;
    }

  output_record ord;
    //ord.rsid=id;
  for(unsigned short int j=0;j<currentPos.size();j++){
      // iterate trough the cursors at every file
      // look for the minimum determined earlier
      // When found, convert p-value
      // add to data structure
      //advance file
    if(currentPos[j]==minimum){
      string tempLine;
        // for each of the cursors that are standing at the minimum
        // push the subset of the vital infos into the "working" arrays
      working_ps.push_back(currentPval[j]);
      working_weights.push_back(weights_p[j]);
      working_rs.push_back(currentPos[j]);
      working_betase.push_back(currentBetaSe[j]);
      working_betas.push_back(currentBeta[j]);
      working_id.push_back(currentRs[j]);
      working_a1.push_back(currentA1[j]);
      working_a2.push_back(currentA2[j]);
      working_af.push_back(currentAF[j]);
        //ord.rsid=currentRs[j];
      working_mask|=(unsigned short int)(pow(2,j));
      if(!getline(*(filestreams[j]), tempLine, '\n')){currentPos[j]="30:1"; info("Read ",counts[j], "lines from file", ifiles[j],".");eofs[j]=true;continue;}

      std::vector <string> line;
        // the following block advances the selected file only (since we have treated the current minimum)
      try{
        counts[j]++;
        line=parse_tab(tempLine, SEP);
        currentPos[j]=line.at(CHR_COL)+":"+line.at(POS_COL);
        currentPval[j]=stold(line[PVAL_COL]);
        currentBetaSe[j]=stold(line[SE_COL]);
        currentBeta[j]=stold(line[BETA_COL]);

        if(ID_COL>-1){currentRs[j]=line[RSID_COL];}else{currentRs[j]=currentPos[j];}
        currentA1[j]=line[A1_COL];
        currentA2[j]=line[A2_COL];
        currentAF[j]=stof(line[AF_COL]);
        if(ID_COL>-1){id=line[RSID_COL];}else{id=line.at(CHR_COL)+":"+line.at(POS_COL);}
        weights_p[j]=SIZE_COL>-1?stoi(line[SIZE_COL]) : sample_sizes[j];
          //cout << "update \t";
          //std::copy(currentPos.begin(), currentPos.end(), std::ostream_iterator<string>(std::cout, " "));
          //cout << "iterator "<<j<<"\n";

      } catch(exception e){
        error("Impossible to parse file ", ifiles[i], 
          "\n\tCheck separator and column numbers for chromosome, position and p-value (", line.size(),")");
      }
    }
  }

    // At this point, working_XX contains the previous minimal positions
    // and current_XX contains the current (next) one

  meta_analyse(working_id, working_rs, working_a1, working_a2, working_ps, working_betas,  working_betase, working_weights, working_mask, working_af, ifiles.size(), correlations, ofs);


  working_mask=0;
  working_weights.clear();
  working_ps.clear();
  working_rs.clear();
  working_betas.clear();
  working_betase.clear();
  working_id.clear();
  working_a1.clear();
  working_a2.clear();
  working_af.clear();

  min= min_element(currentPos.begin(),currentPos.end(),rsid_comparator);
  if(rsid_comparator(*min, minimum)){error("At least one of the input files is unsorted: ", *min, "  <  ", minimum);}
  minimum=*min;                                                                                                                                                                                                                                                      
  if(minimum=="30:1"){break;}
  correlations.update_minima(minimum);

}
info("All done.");
info("METACARPA ( μ - 🐟 ) swimming away.");
info("Goodbye.");

return 0;

}



// DEFUNCT CODE
// OLD BODY OF META_ANALYSE (WITH COMMENTS)


  //   ord.chrpos=working_rs[0];
  //   ord.rsid=working_id[0];
  //   ord.a1=working_a1[0];
  //   ord.a2=working_a2[0];

  //   unsigned int j=0;
  //   for(unsigned short int i=0; i<ifiles.size();i++){
  //     if((working_mask & (1<<i)) == 0){
  //       ord.status.push_back('?');
  //     }else{
  //     //if((mask & mask-1)==0){error("mask is a power of two: "+to_string(mask)+" whereas i="+to_string(i)+" and maskcheck is "+to_string(1<<i));}
  //       char effect_dir=working_betas[j] > 0? '+' : '-';
  //       j++;
  //       ord.status.push_back(effect_dir);
  //     }
  //   }

  //   if(bit_count(working_mask)>1){
  //     // compute weights
  //     int sum_ss=somme(working_weights);
  //     i=0;
  //     //for(long double l : working_weights){working_weights[i]/=sum_ss;i++;}

  //     for(long double l : working_weights){working_weights[i]=sqrt(working_weights[i]);working_weights[i]/=sqrt(sum_ss);i++;}


  //     // compute transforms and sum thereof
  //       i=0;
  //     vector<cpp_dec_float_1000> zs(working_weights.size());
  //     vector<cpp_dec_float_1000> zs_fess(working_weights.size());
  //     for(long double l : working_ps){
  //       zs[i]=z_transform_fess(working_ps[i], working_betas[i],working_rs[i]);
  //       zs_fess[i]=z_transform_fess(working_ps[i], working_betas[i],working_rs[i]);
  //       i++;
  //     }
  //     ord.z=somme(produit(zs, working_weights));


  //     ord.z_fess=somme(produit(zs_fess, working_weights));
  //     // compute p-value SD
  //     // the first product in the line below sums to 1 by definition
  //     //ord.zse=sqrt(somme(produit(working_weights, working_weights))+produit_reciproque_asymetrique(correlations.getmat(working_mask), working_weights));
  //     ord.zse=sqrt(1+produit_reciproque_asymetrique(correlations.getmat(working_mask), working_weights));


  //     //compute matrix for beta weights
  //     //info("Working mask ", working_mask);
  //     ublas::matrix<long double> c=correlations.getmat(working_mask);
  //     for(i=0;i<c.size1();i++){
  //       for(unsigned short int j=0;j<c.size1();j++){
  //         if(i>j){c(i,j)=c(j,i);}else if (i==j){c(i,j)=1;}
  //         c(i,j)=c(i,j)*working_betase[i]*working_betase[j];
  //       }
  //     }
  //     // calculate weights
  //     //    * sum of columns/total sum of matrix
  //     ublas::matrix<long double> inverse(c.size1(),c.size1());

  //     if(!InvertMatrix(c,inverse)){correlations.print(c);error("Could not invert variance/covariance matrix.");}

  //     working_weights=colsum(inverse);
  //     long double totsum=somme(working_weights);
  //     for(i=0;i<working_weights.size();i++) working_weights[i]/=totsum;

  //     // calculate betase
  //       std::vector<long double> working_var=produit(working_betase, working_betase);
  //     ord.betase=sqrt(somme(produit(produit(working_weights, working_weights), working_var)) + produit_reciproque_asymetrique(c, working_weights));

  //     // meta-beta
  //     ord.beta=somme(produit(working_weights, working_betas));
  //         // compute meta-analysis p-value
  //     boost::math::normal_distribution<long double>  correctednormal(0.0, ord.zse);
  //     boost::math::normal_distribution<cpp_dec_float_1000>  correctednormal_p(0.0, ord.zse);
  //     boost::math::normal_distribution<long double>  wald(0.0, 1);
  //     boost::math::normal_distribution<cpp_dec_float_1000>  wald_p(0.0, 1);

  //     try{
  //      ord.p=2*(cdf(correctednormal, -1*abs(static_cast<long double>(ord.z))));
  //      //ord.p_uncorrected=2*(1-cdf(wald, abs(static_cast<long double>(ord.z))));
  //       //ord.p=1-cdf(correctednormal, -1*abs(static_cast<long double>(ord.z)));
  //       // ord.p_uncorrected=1-cdf(wald, -1*abs(static_cast<long double>(ord.z)));
  //      ord.p_fess=2*cdf(wald, -1*abs(static_cast<long double>(ord.z_fess)));
  //      //ord.p=1-cdf(correctednormal, static_cast<long double>(ord.z));
  //       // ord.p_uncorrected=1-cdf(wald, static_cast<long double>(ord.z));
  //      ord.p_wald=2*cdf(wald, -1*abs(ord.beta/ord.betase));
  //      if(ord.p==0){

  //         //info("Warning, position "+minimum+" meta-analyses below long double precision. Analysing with multiprecision...");
  //          //ord.p =1-cdf(correctednormal_p, -1*abs(ord.z));
  //        //ord.p =1-cdf(correctednormal_p, ord.z);
  //         // ord.p_uncorrected=1-cdf(wald_p, abs(ord.z));
  //        ord.p_wald=2*cdf(wald_p, -1*abs(ord.beta/ord.betase));
  //        ord.p=2*(cdf(correctednormal_p, -1*abs(ord.z)));
  //        ord.p_fess=2*cdf(wald_p, -1*abs(ord.z_fess));
  //        //ord.p_uncorrected=2*(1-cdf(correctednormal_p, abs(ord.z)));

  //      }
  //    }catch(exception e){
  //       // ord.p=1-cdf(correctednormal_p, -1*abs(ord.z));
  //     //ord.p=1-cdf(correctednormal_p, ord.z);
  //       // ord.p_uncorrected=1-cdf(wald_p, -1*abs(ord.z));
  //     ord.p=2*(cdf(correctednormal, -1*abs(ord.z)));
  //     //ord.p_uncorrected=2*(1-cdf(wald_p, abs(ord.z)));
  //     ord.p_wald=2*cdf(wald_p, -1*abs(ord.beta/ord.betase));
  //     ord.p_fess=2*cdf(wald_p, -1*abs(ord.z_fess));
  //   }

  // }else{
  //     // SNP present in only 1 study, skip all calculations
  //   ord.rsid=working_id[0];
  //   ord.beta=working_betas[0];
  //   ord.betase=working_betase[0];
  //   ord.p=-1;
  //   ord.p_wald=-1;
  //   ord.p_fess=-1;
  //   ord.z=-1;
  //   ord.zse=-1;
  //   ord.p_uncorrected=-1;
  //   ord.a1=working_a1[0];
  //   ord.a2=working_a2[0];
  // }
  // ord.print(ofs);
  /*
  map <string, position_info, cmp_rsid> thisstudy;
  bool firstline=true;
  info("Calculating z-scores and binomial transformation for file "+filename+"...");
  string tempLine;
  boost::math::normal_distribution<long double>  standardnormal;
  boost::math::normal_distribution<cpp_dec_float_1000>  standardnormal_p;
  means[i]=0;
  while ( getline(inputFile, tempLine, '\n') ) {
    if(firstline==true){firstline=false;continue;}
    vector <string> line;
    stringstream ss(tempLine);
    string temp;
    while (getline(ss, temp, '\t')) {  
      line.push_back(temp);
    }
    string chrpos=line[0]+":"+line[2];
    position_info thisposition;
      //thisposition.beta=stold(line[7]);
      //means[i]=means[i]+thisposition.beta;
      //meancount[i]++;
      //thisposition.sdbeta=stold(line[8]);
    thisposition.pvalue=stold(line[13]);
    thisposition.rsid=line[1];
    try
    {
      thisposition.zf=(cpp_dec_float_1000)(quantile(standardnormal, 1-thisposition.pvalue/2)*sign(thisposition.beta));
      thisposition.bpval=(thisposition.zf<=0);
      thisposition.need_precision=false;
    }catch(exception &ex){
      thisposition.need_precision=true;
      info("\tPosition "+chrpos+" could not be transformed with normal machine precision. Using arbitrary precision (p="+line[13]+")");
      cpp_dec_float_1000 p1=1-(cpp_dec_float_1000)thisposition.pvalue/2;
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



  // All studies have been read. Do we have positions only present in one study?
// typedef map <string, unsigned short int, cmp_rsid>::iterator it_type;
// for(it_type it =where_present.begin(); it != where_present.end();){
//   unsigned short int mask=it->second;
//   if((mask & mask-1)==0){removed++;}
//   ++it;
// }
// info(to_string(removed)+" positions were only genotyped in one study. Sadly, they will be lost.");


// Time to compute relatedness matrix.
// This calculates the correlations using the binomially tansformed betas.
// It is the main innovation introduced by Province and Borecki
info("Calculating p-value correlation matrix...");
vector<vector<long double>> cormat(ifiles.size(), vector<long double>(ifiles.size(), 0));
i=0;
vector<long double> sds(ifiles.size(), 0);
vector<bool> done(ifiles.size(), false);
vector<int> counts(ifiles.size(), 0);
for_each(ifiles.begin(), ifiles.end(), [&](string s)-> void {
  unsigned short int j=i+1;
  unsigned short int expj=1<<j;
  unsigned short int expi=1<<i;
  for_each(next(ifiles.begin(),i+1), ifiles.end(), [&](string t)-> void {
        // get all positions that have at least these two studies in common
    info("between file "+s+" and file "+t);
    if (i==j){cormat[i][j]=1;}else{
      typedef map <string, unsigned short int, cmp_rsid>::iterator it_type;
      vector<bool> a;
      vector<bool> b;
      for(it_type it =where_present.begin(); it != where_present.end();){
        unsigned short int mask=it->second;
        string chrpos=it->first;
        // Here we calculate SDs for entire studies...
        //      info("Treating position "+chrpos+" with mask "+to_string(mask));
        if(!done[j] && (mask & expj)==expj){sds[j]=sds[j]+(studydata[j][chrpos].beta-means[j])*(studydata[j][chrpos].beta-means[j]);counts[j]++;}
        if(!done[i] && (mask & expi)==expi){sds[i]=sds[i]+(studydata[i][chrpos].beta-means[i])*(studydata[i][chrpos].beta-means[i]);counts[i]++;}
        if((mask & expj)==expj && (mask & expi)==expi){
        //            info("Position "+chrpos+" is present in both studies.");
          a.push_back(studydata[i][chrpos].bpval);
          b.push_back(studydata[j][chrpos].bpval);

        }
        ++it;
      }
      done[j]=true;
      done[i]=true;

      cormat[i][j]=tetrachoric(a,b);

    }
    info("\t...is "+to_string(cormat[i][j]));
    j++;
  });
i++;
});
studies_correlation sc=studies_correlation();
sc.add_mask((unsigned short int)3);
bool v[]={true, true};
sc.update_mask((unsigned short int)3, v);
v[1]=false;
sc.update_mask((unsigned short int)3, v);
v[0]=false;
sc.update_mask((unsigned short int)3, v);
v[0]=true;
sc.update_mask((unsigned short int)3, v);
v[1]=true;
sc.update_mask((unsigned short int)3, v);
v[0]=false;
sc.update_mask((unsigned short int)3, v);
v[1]=false;
sc.update_mask((unsigned short int)3, v);
sc.update_mask((unsigned short int)3, v);
v[0]=true;
v[1]=true;
sc.update_mask((unsigned short int)3, v);
v[0]=false;
sc.update_mask((unsigned short int)3, v);
v[1]=true;
v[0]=true;
sc.update_mask((unsigned short int)3, v);
v[0]=false;
sc.update_mask((unsigned short int)3, v);
v[1]=false;
sc.update_mask((unsigned short int)3, v);
sc.update_mask((unsigned short int)3, v);




sc.compute_correlations();

sc.write("test_serial");

studies_correlation scr=studies_correlation();
scr.read("test_serial");
info(scr.getmat((unsigned short int)3));
*/





















