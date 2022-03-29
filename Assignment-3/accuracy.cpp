#include <bits/stdc++.h>
using namespace std;

#define oper remove_if(tokenized.begin(), tokenized.end(),[](string const& s) { return s.size() == 0;})


vector<string>tokenize(const string str,const regex re){
    int c=count(str.begin(),str.end(),',');    
	sregex_token_iterator it{ str.begin(), 
							str.end(), re, -1 }; 
	vector<string> tokenized{ it, {} }; 

	
	tokenized.erase(oper,tokenized.end());
	return tokenized; 
}

void fileread(string s,vector<set<int>> &V){    
    const regex re(R"([\s|,]+)");    
    ifstream in(s);
    string line;
    while(getline(in,line)){
        vector<string> temp=tokenize(line,re);
        set<int>vec;
        for(auto x:temp){
            vec.insert(stoi(x));
        }
        V.push_back(vec);     

    return;
    }

}

int main(int argc, char* argv[]){

    string gt=argv[1];
    string output_file=argv[2];
    int k=stoi(argv[3]);

    vector<set<int>>Ground_truth;
    vector<set<int>>Output;
    fileread(gt,Ground_truth);
    fileread(output_file,Output);

    double avg_prec=0.0;
    double avg_recall=0.0;

    for(int i=0;i<Ground_truth.size();i++){
        vector<int>common_data;
        set<int>s1=Ground_truth[i];
        set<int>s2=Output[i];
        set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),back_inserter(common_data));
        int intersect_size=common_data.size();
        double t1;
        double t2;
        t1=(double)(intersect_size)/(double) k;
        t2=(double)(intersect_size)/(double)(Ground_truth[i].size());
        avg_prec+=t1;
        avg_recall+=t2;
    


}

avg_prec=avg_prec/(double)Ground_truth.size();
avg_recall=avg_recall/(double)Ground_truth.size();

cout<<"Precision: "<<avg_prec<<"\n";
cout<<"Recall: "<<avg_recall<<"\n";

}

