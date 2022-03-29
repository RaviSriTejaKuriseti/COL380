#include <bits/stdc++.h>
using namespace std;

int bytes_to_int(char* arr){
union{
    int f;
    char buffer[4];
}u;

memcpy(&u.buffer,arr,4);
return u.f;
}



int main(int argc, char* argv[]){


ifstream fin;
    fin.open("helper_n.txt");
    int num_usr;
    int k;
    string out_file;
    fin>>num_usr>>k>>out_file;
    fin.close();


//string out_file="user_pred.txt";
const char *fname="output.dat";
FILE *fp=fopen(fname,"rb");
ofstream fout;
fout.open(out_file);
char bytes[4];
int val=-1;
for(int i=0;i<num_usr;i++){
    for(int j=0;j<k-1;j++){
        fread(bytes,4,1,fp);
        val=bytes_to_int(bytes);
        fout<<val<<" ";

    }
    fread(bytes,4,1,fp);
    val=bytes_to_int(bytes);
    fout<<val<<"\n";

}
fout.close();
fclose(fp);    




}