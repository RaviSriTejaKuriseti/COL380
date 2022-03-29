//use this to read binary files from path and execute code


#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>


using namespace std;




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


void read_vect_file(string vec_file,int rows,int cols,vector<vector<float>> &V){
    const char *fname=vec_file.c_str();

    FILE *fp=fopen(fname,"rb");
    char bytes[4];
    float val;

    for(int i=0;i<rows;i++){
        vector<float>vf;
        for(int j=0;j<cols;j++){
            fread(bytes,4,1,fp);
            val=bytes_to_float(bytes);
            vf.push_back(val);
        }
        V.push_back(vf);  
    }
    fclose(fp);
}


void read_int_file(string vec_file,int l,vector<int> &V){
    const char *fname=vec_file.c_str();

    FILE *fp=fopen(fname,"rb");
    char bytes[4];
    int val;

    for(int i=0;i<l;i++){
        fread(bytes,4,1,fp);
        val=bytes_to_int(bytes);
        V.push_back(val);
    }  
    
    fclose(fp);
}





int read_user_file(string s,int cols,vector<vector<float>> &U){
    int rows=0;
    float a;
    ifstream fin(s);
    if(fin.is_open()){
        while(fin){
            vector<float>vec;
            fin>>a;
            if(fin.eof()){
                break;
            }
            vec.push_back(a);           
            for(int i=1;i<cols;i++){
                fin>>a;
                vec.push_back(a);
            }
           U.push_back(vec);           
           rows+=1;
        }
    }
    fin.close();
    cout<<rows<<" "<<cols<<"\n";
    return rows;
   

}


void write_to_file(string in_file,string out_file,int usr_num,int k){
    const char *fname=in_file.c_str();
    FILE *fp=fopen(fname,"rb");
    ofstream fout;
    fout.open(out_file);
    char bytes[4];
    int val=-1;
    for(int i=0;i<usr_num;i++){
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
    return;
   

}




struct lowers{ 
    bool operator()(const pair<int,float>& a,const  pair<int,float>& b) const{
        if(a.second!=b.second){
            return a.second<b.second;
        }
        return a.first<b.first;  
    } 
};


float cosine_dist(vector<float> &vec1,vector<float> &vec2){
    if(vec1.size()!=vec2.size()){
        return -2.00;
    }
    float s1=0;
    float s2=0;
    float s3=0;
    for(int i=0;i<vec1.size();i++){
        float t1=vec1[i];
        float t2=vec2[i];
        s1+=(t1*t1);
        s2+=(t2*t2);
        s3+=(t1*t2);
    }
    float pro1=sqrt(s1*s2);
    float ans=1-s3/pro1;
    return ans;

}


float get_max(vector<pair<int,float>> &vec){
    float x1=vec.front().second;
    for(auto x:vec){
        if(x.second>x1){
            x1=x.second;
        }
    }
    return x1;
}






vector<pair<int,float>> SearchLayer(int k,vector<float>&q, vector<pair<int,float>> &candidates,
 vector<int>&indptr,vector<int>&index,vector<int>&level_offset,int lc, map<int,bool> &visited,vector<vector<float>> &vect){


vector<pair<int,float>>top_k;
top_k.insert(top_k.begin(),candidates.begin(),candidates.end());

while(candidates.size()>0){
    pop_heap(candidates.begin(), candidates.end(),lowers());
    pair<int,float> pep=candidates.back();
    int ep=pep.first;
    candidates.pop_back();
    int start=indptr[ep]+level_offset[lc];
    int end=indptr[ep]+level_offset[lc+1];

    for(int j=start;j<end;j++){
        auto itr=visited.find(index[j]);
        if(itr!=visited.end() || index[j]==-1){
            continue;
        }
        visited.insert({index[j],true});
        float _dist=cosine_dist(q,vect[index[j]]);

        float max_dist=get_max(top_k);
        if(_dist>max_dist && top_k.size()>=k){
            continue;
        }

        top_k.push_back(make_pair(index[j],_dist));

        sort(top_k.begin(),top_k.end(),lowers());
        if(top_k.size()>=k){
            top_k.resize(k);
        }
        make_heap(top_k.begin(),top_k.end(),lowers());

        candidates.push_back(make_pair(index[j],_dist));
        push_heap(candidates.begin(),candidates.end(),lowers()); 
       
    }
    
    
}


return top_k;
}



vector<pair<int,float>> QueryHNSW(int k,vector<float>&q,int ep,vector<int>&indptr,vector<int>&index,vector<int>&level_offset,int max_level,vector<vector<float>> &vect){

    vector<pair<int,float>>top_k;
    top_k.push_back(make_pair(ep,cosine_dist(q,vect[ep])));
    push_heap(top_k.begin(),top_k.end(),lowers());
    map<int,bool>visited;
    visited.insert({ep,true});
    for(int lvl=max_level;lvl>=0;lvl--){
        top_k=SearchLayer(k,q,top_k,indptr,index,level_offset,lvl,visited,vect);
    }
    sort(top_k.begin(),top_k.end(),lowers()); 
    return top_k;
}







void get_val(vector<int> &A,vector<pair<int,float>> &V){
    for(auto x:V){
        A.push_back(x.first);
    }
}






int main(int argc, char* argv[]){
    //assert(argc > 4);
    string out_path = argv[1]; 
    //read from here ep,max_level,dims(in one file),indptr.dat,index.dat,level_offset.dat,vect.dat
    int k=stoi(argv[2]);
    //value of k
    string user_file = argv[3];
    //read user.txt
    string output_file=argv[4];
    //write predictions.txt

    ifstream fin;
    fin.open(out_path+"/"+"helper.txt");
    int ep;
    int max_level;
    int rows;
    int cols;
    int si;
    int sptr;
    int slo;
    fin>>max_level>>ep>>rows>>cols>>si>>sptr>>slo;
    fin.close();
 
   

    vector<vector<float>>U; //user-vector.
    int num_usr=read_user_file(user_file,cols,U);  //number of users


    ofstream ff;
    ff.open("helper_n.txt");
    ff<<num_usr<<" "<<k<<" "<<output_file;
    ff.close();
 

    vector<vector<float>>vec;
    read_vect_file(out_path+"/"+"vect.dat",rows,cols,vec);

    vector<int>indptr;
    vector<int>level_offset;
    vector<int>index;

    read_int_file(out_path+"/"+"indptr.dat",sptr,indptr);
    read_int_file(out_path+"/"+"level_offset.dat",slo,level_offset);
    read_int_file(out_path+"/"+"index.dat",si,index);

    


    int rank,size;  //size is the number of threads
    
    //Starting MPI pipeline
    MPI_Init(NULL, NULL);
    
    // // Extracting Rank and Processor Count
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int thread_num=2;
    bool flag=true;

    if(rank==0){
        ofstream fw1("output.dat",ios::binary|ios::out);
        fw1.seekp(4*(num_usr)*k-1);
        fw1.write("",1);
        fw1.close();
        for(int I=1;I<size;I++){                
            MPI_Send(&flag,1,MPI_INT,I,0,MPI_COMM_WORLD);
        }
             
    }


    if(rank!=0){
        MPI_Recv(&flag,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }


    ofstream fw("output.dat",ios::binary|ios::out|ios::in);
    cout<<rank<<" "<<"\n";


    // #pragma omp
    int n2=((rank+1)*num_usr)/size;
    int n1=(rank*num_usr)/size;

    int* arr=new int[(n2-n1)*k];

    #pragma omp parallel
    {

        #pragma omp for

            for(int i=(rank*num_usr)/size;i<((rank+1)*num_usr)/size;i++){

             //int tid=omp_get_thread_num();

               //cout<<rank<<" "<<tid<<"\n";
               //vector<pair<int,float>>top_k;
               vector<pair<int,float>>ans;
               vector<int>V;
               ans=QueryHNSW(k,U[i],ep,indptr,index,level_offset,max_level,vec);
               get_val(V,ans);
               for(int j=0;j<k;j++){
                   arr[(i-n1)*k+j]=V[j];
               }
               
               //Here write the ans to output_file

            }
        


    }
    int temp=0;
    fw.seekp(k*4*n1,ios::beg);
    char bytes[4];
    for(int i=n1;i<n2;i++){
        for(int j=0;j<k;j++){
            int x=arr[temp];
            int_to_bytes(x,bytes);
            fw.write(bytes,4);
            temp+=1;
        }
        
    }
    delete arr;

    fw.close();


    
    MPI_Finalize();
   
}
