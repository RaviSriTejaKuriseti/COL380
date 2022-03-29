#include <bits/stdc++.h>
#define oper remove_if(tokenized.begin(), tokenized.end(),[](string const& s) { return s.size() == 0;})
using namespace std;

void float_to_bytes(float val,char* arr){
union{
    float f;
    char buffer[4];
}  u;

u.f=val;
memcpy(arr,&u.buffer,4);

}


float bytes_to_float(char* arr){
union{
    float f;
    char buffer[4];
}  u;

memcpy(&u.buffer,arr,4);
return u.f;
}


void int_to_bytes(int val,char* arr){
    union{
    int f;
    char buffer[4];
}u;

u.f=val;
memcpy(arr,&u.buffer,4);


}

int bytes_to_int(char* arr){
union{
    int f;
    char buffer[4];
}u;

memcpy(&u.buffer,arr,4);
return u.f;
}





vector<string>tokenize(const string str,const regex re){
     int c=count(str.begin(),str.end(),',');    
	sregex_token_iterator it{ str.begin(), 
							str.end(), re, -1 }; 
	vector<string> tokenized{ it, {} }; 

	
	tokenized.erase(oper,tokenized.end()); 
	return tokenized; 
}


pair<int,int> vect_file_io(string s,string out_file){  
    const regex re(R"([\s|,]+)");
    int cols=0; 
    int rows=0; 
    ifstream in(s);
    ofstream fout(out_file,ios::binary|ios::out);
    string line;
    float a;
    char bytes[4];
    getline(in,line);
    in.close();
    vector<string> temp=tokenize(line,re);
    cols=temp.size(); 
    in.close();   

    ifstream fin(s);
    if(fin.is_open()){
        while(fin){
            fin>>a;
            if(fin.eof()){
                break;
            }
            float_to_bytes(a,bytes);
            fout.write(bytes,4);
            for(int i=1;i<cols;i++){
                fin>>a;
                float_to_bytes(a,bytes);
                fout.write(bytes,4);
            }
           
           rows+=1;
        }
    }
    fin.close();
    fout.close();


    pair<int,int>P;
    P.first=rows;
    P.second=cols;
    return P;
   

}


int write_file_to_bin(string in_path,string out_path,string fname){

    int l;
    char bytes[4];
    int ct=0;

    ifstream fin;
    string s;
    string s1=fname.substr(0,fname.length()-3);
    s=s1+"dat";
    ofstream fout(out_path+"/"+s ,ios::binary|ios::out);
    fin.open(in_path+"/"+fname);
     if(fin.is_open()){
        while(fin){
            fin>>l;
            if(fin.eof()){
                break;
            }
            int_to_bytes(l,bytes);
            fout.write(bytes,4);
            ct+=1;
        }
    }
    fin.close();
    fout.close();
    return ct;

}





void fileread(string infile_path,vector<vector<float>> &vec){     
    ifstream in(infile_path);
    string line;
    int count=0;
    while(getline(in,line)){
        cout<<line<<"\n";
        vector<string> temp;
    }
   
    return;

}





int main(int argc, char* argv[]){
    string in_path =argv[1];
    string out_path=argv[2];

    //write_file_to_bin(in_path,out_path,"level.txt");
    int si=write_file_to_bin(in_path,out_path,"index.txt");
    int sptr=write_file_to_bin(in_path,out_path,"indptr.txt");
    int slo=write_file_to_bin(in_path,out_path,"level_offset.txt");

    ifstream fin;
    ofstream fout;
    string vect_file=in_path+"/"+"vect.txt";
    string vect_out_file=out_path+"/"+"vect.dat";

    char bytes[4];
    pair<int,int> P=vect_file_io(vect_file,vect_out_file);


    //cout<<"dat file writing done"


    
    fin.open(in_path+"/"+"max_level.txt");
    int lvl;
    fin>>lvl;
    fin.close();
    fin.open(in_path+"/"+"ep.txt");
    int ep;
    fin>>ep;
    fin.close();
    fout.open(out_path+"/"+"helper.txt");
    fout<<lvl<<" "<<ep<<" "<<P.first<<" "<<P.second<<" "<<si<<" "<<sptr<<" "<<slo<<"\n";
    fout.close();
    


}