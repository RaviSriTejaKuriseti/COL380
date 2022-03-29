#include <bits/stdc++.h>
#include <string>
#include <mpi.h>
#include <assert.h>
#include "randomizer.hpp"

using namespace std;


//Notice how the randomizer is being used in this dummy function
void print_random(int tid, int num_nodes, Randomizer r){
    int this_id = tid;
    int num_steps = 10;
    int num_child = 20;

    std::string s = "Thread " + std::to_string(this_id) + "\n";
    std::cout << s;

    for(int i=0;i<num_nodes;i++){
        //Processing one node
        for(int j=0; j<num_steps; j++){
            if(num_child > 0){
                //Called only once in each step of random walk using the original node id 
                //for which we are calculating the recommendations
                int next_step = r.get_random_value(i);
                //Random number indicates restart
                if(next_step<0){
                    std::cout << "Restart \n";
                }else{
                    //Deciding next step based on the output of randomizer which was already called
                    int child = next_step % 20; //20 is the number of child of the current node
                    std::string s = "Thread " + std::to_string(this_id) + " rand " + std::to_string(child) + "\n";
                    std::cout << s;
                }
            }else{
                std::cout << "Restart \n";
            }
        }
    }
}


bool comp(pair<int,int> &P1,pair<int,int> &P2){
    if(P1.second!=P2.second){
        return P1.second>P2.second;
    }
    return P1.first<P2.first;
    
    
    

}





int main(int argc, char* argv[]){
    assert(argc > 8);
    std::string graph_file = argv[1];
    int num_nodes = std::stoi(argv[2]);
    int num_edges = std::stoi(argv[3]);
    float restart_prob = std::stof(argv[4]);
    int num_steps = std::stoi(argv[5]);
    int num_walks = std::stoi(argv[6]);
    int num_rec = std::stoi(argv[7]);
    int seed = std::stoi(argv[8]);

    
    //Only one randomizer object should be used per MPI rank, and all should have same seed
    Randomizer random_generator(seed, num_nodes, restart_prob);
    int rank, size;


    auto begin = std::chrono::high_resolution_clock::now();

    vector<int>out_degrees(num_nodes,0);
    vector<vector<int>>G;
    for(int i=0;i<num_nodes;i++){
        vector<int>vec;
        G.push_back(vec);
    }

    vector<map<int,int>>inf_score;
    for(int i=0;i<num_nodes;i++){
        map<int,int>M;
        inf_score.push_back(M);
    }

    const char *fname=graph_file.c_str();

    FILE *fp=fopen(fname,"rb");
    unsigned char bytes[4];
    

    for(int i=0;i<num_edges;i++){        
        int s1=0;
        int s2=0;
        fread(bytes,4,1,fp);
        s1+=bytes[3]|(bytes[2]<<8)|(bytes[1])<<16|(bytes[0]<<24);
        fread(bytes,4,1,fp);
        s2+=bytes[3]|(bytes[2]<<8)|(bytes[1])<<16|(bytes[0]<<24);    
        out_degrees[s1]+=1;
        G[s1].push_back(s2);
    }
    fclose(fp);

       

    

      

    
    //Starting MPI pipeline
    MPI_Init(NULL, NULL);
    
    // // Extracting Rank and Processor Count
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    

   


    int nodes_per_process=num_nodes/size;
    int rem=num_nodes % size;
    int start_ind=0;
    int end_ind=0;
    vector<int>Start(size,0);
    vector<int>End(size,0);
   

    if(rem==0){
       start_ind=rank*nodes_per_process;
       end_ind=start_ind+nodes_per_process;
    }
    else{
        start_ind=rank*(nodes_per_process+1);
        end_ind=min(start_ind+nodes_per_process+1,num_nodes);

    }   
    

    Start[rank]=start_ind;
    End[rank]=end_ind;
    bool flag=true; 

    if(rank==0){
        ofstream fw1("output.dat",ios::binary|ios::out);
        fw1.seekp((4+8*num_rec)*num_nodes-1);
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
    fw.seekp(start_ind*(2*num_rec+1)*4 + ios::beg, ios::beg);
    for(int i=start_ind;i<end_ind;i++){
        for(int j=0;j<G[i].size();j++){            
            for(int k=0;k<num_walks;k++){
                int curr_node=G[i][j];
                int v=0;
                while(v<num_steps){
                    int l=G[curr_node].size();
                    if(l==0){                        
                        curr_node=G[i][j];
                    }
                    else{
                        int nxt_child=random_generator.get_random_value(i);
                        if(nxt_child==-1){                        
                            curr_node=G[i][j];
                        }
                        else{
                            curr_node=G[curr_node][nxt_child%l];
                           
                        }

                        
                    }
                    auto itr=inf_score[i].find(curr_node);
                    if(itr==inf_score[i].end()){
                        inf_score[i].insert(pair<int,int>(curr_node,1));
                    }
                    else{
                        itr->second+=1;
                    }
                    v+=1;

        

                }
                

            }
            
        }




        vector<pair<int,int>>PV;

         for(int u=0;u<G[i].size();u++){
            inf_score[i].erase(G[i][u]);
        }
        inf_score[i].erase(i);
        for(auto x:inf_score[i]){
            PV.push_back(x);
        }  
   
        sort(PV.begin(),PV.end(),comp);
        unsigned char bytes[4];
        int temp=out_degrees[i];

        bytes[0]=(temp>>24)&255;
        bytes[1]=(temp>>16)&255;
        bytes[2]=(temp>>8)&255;
        bytes[3]=temp&255;
        fw.write((char*) bytes,4);

        int l1=PV.size();

        for(int u=0;u<min(num_rec,l1);u++){
            pair<int,int>P=PV[u];            
            temp=P.first;
            bytes[0]=(temp>>24)&255;
            bytes[1]=(temp>>16)&255;
            bytes[2]=(temp>>8)&255;
            bytes[3]=temp&255;
            fw.write((char*) bytes,4);
            temp=P.second;
            bytes[0]=(temp>>24)&255;
            bytes[1]=(temp>>16)&255;
            bytes[2]=(temp>>8)&255;
            bytes[3]=temp&255;
            fw.write((char*) bytes,4);
            
        }
                
        
        for(int u=min(num_rec,l1);u<num_rec;u++){
            bytes[0]='N';
            bytes[1]='U';
            bytes[2]='L';
            bytes[3]='L';
            fw.write((char*) bytes,4);
            fw.write((char*) bytes,4);
        }

    }
    fw.close();  
    
    //cout<<rank<<"\n";
            
        
     
    MPI_Finalize();
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    
    double duration = (1e-6 * (std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin)).count());
    if(rank==0){
        std::cout << "Time taken for process: " << duration << "ms" << "\n";
    }
    
}