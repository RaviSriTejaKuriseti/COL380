// %%cuda --name hello.cu

#include <bits/stdc++.h>
using namespace std;



__device__ float sine(int x){
    float INV_ROOT_TWO= (1.0f)/sqrt(2.0f);
    if(x==0){
        return 0;
    }
    else if(x==45){
        return INV_ROOT_TWO;
    }
    else{
        return  -1*INV_ROOT_TWO;
    }
}


__device__ float cosine(int x){
    float INV_ROOT_TWO= (1.0f)/sqrt(2.0f);
    if(x==0){
        return 1;
    }
    else{
        return INV_ROOT_TWO;
    }
}



__device__ void get_new_coord(float* point,float* pivot,float* new_pt,int theta){  //[0] is x-coord and [1] is y-coord
    //left corner is origin rightwards +ve x-axis and top +ve y-axis 
    new_pt[0]=pivot[0]+(point[0]-pivot[0])*cosine(theta)-(point[1]-pivot[1])*sine(theta);
    new_pt[1]=pivot[1]+(point[0]-pivot[0])*sine(theta)+(point[1]-pivot[1])*cosine(theta);
    return;    


}




__device__ void get_bound_pos(int r,int c,int theta,int *bdrs){

    float m1=(float) c;
    float n1=(float) r;



    float a1[2]={m1,0.0f}; //bottom right
    float a2[2]={m1,n1}; //top right
    float a3[2]={0.0f,n1};  //top left
    float P[2]={0.0f,0.0f};

    float p1[2];
    float p2[2];
    float p3[3];


    get_new_coord(a1,P,p1,theta); //new pos of bottom right
    get_new_coord(a2,P,p2,theta); //new pos of top right
    get_new_coord(a3,P,p3,theta); //new pos of top left

    int lx,rx,ty,by;  //left-x,right-x,top-y,bottom-y

    if(theta==0){
        lx=0;
        rx=c;
        ty=r;
        by=0;

    }
    else if(theta==45){
        lx=ceil(p3[0]);
        rx=floor(p1[0]);
        ty=floor(p2[1]);
        by=0;


    }
    else if(theta==-45){
        lx=0;
        rx=floor(p2[0]);
        ty=floor(p3[1]);
        by=ceil(p1[1]);


    }

    bdrs[0]=lx;
    bdrs[1]=rx;
    bdrs[2]=ty;
    bdrs[3]=by;



    return;
  

}


__device__ bool check_in_limits(float* Grey_img,int* bounds,int d_m,int d_n,int i,int j,float th2_val,float avg_val,float* d_A){

    int l1=d_m;
    int l2=d_n;

    int nl1=bounds[1]-bounds[0]+1;
    int nl2=bounds[2]-bounds[3]+1;
    float val=0.0f;

    if(j+bounds[0]<0 || j+bounds[1]>=l2 || i+bounds[2]>=l1 || i+bounds[3]<0){
        return false;
    }

    for(int u=i+bounds[3];u<=i+bounds[2];u++){
        for(int v=j+bounds[0];v<=j+bounds[1];v++){
           
            val+=Grey_img[u*d_n+v];

        }

        
    }
    val=(val)/(float)(nl1*nl2);
  
    if(abs(val-avg_val)<th2_val){
        d_A[i*d_n+j]=val-avg_val;
        //cout<<i<<" "<<j<<" "<<val<<"\n";
        return true;
    }
    d_A[i*d_n+j]=-1.0f;
    return false;


}

__device__ bool interpolate(float* point,float* val,int* D,int d_m,int d_n){ //vector<float> &point,vector<vector<vector<int>>> &D){
    int x1=floor(point[0]);
    int y1=floor(point[1]);
    int x2=1+floor(point[0]);
    int y2=1+floor(point[1]);
    
    

     if(x1<0 || y1<0 || x2>=d_n || y2>=d_m){
        return false;

     }

   

    //z00*(1-x)*(1-y) + z10*x*(1-y) + z01*(1-x)*y + z11*x*y

    float z00;
    float z10;
    float z01;
    float z11;



      
    
    for(int u=0;u<3;u++){
       
        z00=(float) D[y1*3*d_n+x1*3+u];
        z01=(float) D[y2*3*d_n+x1*3+u];
        z10=(float) D[y1*3*d_n+x2*3+u];
        z11=(float) D[y2*3*d_n+x2*3+u];

        val[u]=z00*(x2-point[0])*(y2-point[1])+z01*(point[1]-y1)*(x2-point[0])+z10*(point[0]-x1)*(y2-point[1])+z11*(point[0]-x1)*(point[1]-y1);
    }



    return true;
}


__device__ void get_best_picks(int* D,int* Q,int i,int j,float th1_val,int theta,int q_m,int q_n,int d_m,int d_n,int top_n,float* d_R){


    float point[2];
    float pivot[2];
    float new_pos[2];
    float interpol_val[4];

    point[0]=0.0f;
    point[1]=0.0f;

    pivot[0]=(float)j;
    pivot[1]=(float)i;
    

    float rmsd=0.0f;
    float temp;
    float f;
    bool pol_flag;

    for(int u=i;u<i+q_m;u++){
        for(int v=j;v<j+q_n;v++){
            point[0]=(float)v;
            point[1]=(float)u;
            get_new_coord(point,pivot,new_pos,theta);
            pol_flag=interpolate(new_pos,interpol_val,D,d_m,d_n);
            if(pol_flag==false){
                d_R[i*d_n+j]=-1.0f;
                return;
            }
            for(int w=0;w<3;w++){
                temp=(Q[(u-i)*3*q_n+(v-j)*3+w]-interpol_val[w]);
                rmsd+=temp*temp;

            }

           

        }


    }

    float n1=(float) (3*q_m*q_n);
    rmsd/=n1;
    f=sqrt(rmsd);
    //cout<<i<<" "<<j<<" "<<theta<<" "<<f<<"\n";
    if(f<th1_val){
        d_R[i*d_n+j]=f;
        //cout<<i<<" "<<j<<" "<<theta<<" "<<"\n";
        return;

    }
    d_R[i*d_n+j]=-1.0f;


   
    return;
   
}
   







__global__ void filtering(float* Grey_img,int theta,float th2_val,float avg_val,int q_m,int q_n,int d_m,int d_n,
int* D,int* Q,float th1_val,int top_n,float* d_A,float* d_R){

    int bounds[4];
    get_bound_pos(q_m,q_n,theta,bounds);
    int i=blockIdx.y;
    int j=blockIdx.x;


    bool filt_flag=check_in_limits(Grey_img,bounds,d_m,d_n,i,j,th2_val,avg_val,d_A);
    if(filt_flag){
       get_best_picks(D,Q,i,j,th1_val,theta,q_m,q_n,d_m,d_n,top_n,d_R);              

    }
    else{
        d_R[i*d_n+j]=-1.0f;
    }
    

    return;


}







int main(int argc,char* argv[]){

    auto beg=std::chrono::high_resolution_clock::now();

    string data_img_path=argv[1];
    string query_img_path=argv[2];
    float th_1=stof(argv[3]);
    float th_2=stof(argv[4]);
    int n=stoi(argv[5]);

    

    int q_m,q_n;
    int d_m,d_n;

    int s=0;
    float avg_val;
    int val;

    
    ifstream fin;
   
    fin.open(data_img_path);
   
    fin>>d_m>>d_n;
    int* D=new int[3*d_m*d_n];
    
    for(int i=0;i<d_m;i++){
        for(int j=0;j<d_n;j++){
            for(int k=0;k<3;k++){
                fin>>val;
                D[(d_m-1-i)*d_n*3+j*3+k]=val;
            }
        }

    }
  
    fin.close();
    fin.open(query_img_path);
    cout<<query_img_path<<"\n";
   
    fin>>q_m>>q_n;
    int* Q=new int[3*q_m*q_n];
      
    for(int i=0;i<q_m;i++){
        for(int j=0;j<q_n;j++){
            for(int k=0;k<3;k++){
                fin>>val;
                Q[(q_m-1-i)*q_n*3+j*3+k]=val;
                s+=val;
            }
        }

    }
    
    fin.close();

   
    avg_val=(float)(s)/(float)(3*q_m*q_n);

    float ans=0.0f;
    float* G=new float[d_m*d_n];

    for(int i=0;i<d_m*d_n;i++){
        ans=(float)(D[3*i]+D[3*i+1]+D[3*i+2])/(3.0f);
        G[i]=ans;
    }
    
    
    cout<<q_m<<" "<<q_n<<" "<<d_m<<" "<<d_n<<"\n";


    int *d_Q;
    int *d_D;
    float *d_G;

    float* d_A;
    float* d_R;

    float* A=new float[d_m*d_n];
    float* R=new float[d_m*d_n];




    cudaMalloc((void**)&d_Q,sizeof(float)*3*q_m*q_n);
    cudaMalloc((void**)&d_D,sizeof(float)*3*d_m*d_n);
    cudaMalloc((void**)&d_G,sizeof(float)*d_m*d_n);


    cudaMalloc((void**)&d_A,sizeof(float)*d_m*d_n);
    cudaMalloc((void**)&d_R,sizeof(float)*d_m*d_n);


    cudaMemcpy(d_D, D,sizeof(float)*3*d_m*d_n,cudaMemcpyHostToDevice);
    cudaMemcpy(d_G, G, sizeof(float)*d_m*d_n, cudaMemcpyHostToDevice);
    cudaMemcpy(d_Q, Q, sizeof(float)*3*q_m*q_n, cudaMemcpyHostToDevice);


    delete D;
    delete G;
    delete Q;



    
    


    priority_queue<pair<float,tuple<int,int,int,float>>>PQ;
    vector<pair<float,tuple<int,int,int,float>>>V;
    vector<int>angles{-45,0,45};

    dim3 grid(d_n,d_m);

    for(auto h:angles){
        filtering<<<grid,1>>>(d_G,h,th_2,avg_val,q_m,q_n,d_m,d_n,d_D,d_Q,th_1,n,d_A,d_R);
        cudaMemcpy(A,d_A, sizeof(float)*d_m*d_n, cudaMemcpyDeviceToHost);
        cudaMemcpy(R,d_R, sizeof(float)*d_m*d_n, cudaMemcpyDeviceToHost);
       

        for(int i=0;i<d_m;i++){
            for(int j=0;j<d_n;j++){
                if(R[i*d_n+j]>=0.0f){
                    pair<float,tuple<int,int,int,float>>P;
                    P.first=R[i*d_n+j];
                    get<0>(P.second)=i;
                    get<1>(P.second)=j;
                    get<2>(P.second)=h;
                    get<3>(P.second)=A[i*d_n+j];
                    PQ.push(P);
                    if(PQ.size()>n){
                        PQ.pop();
                    }
                }
            }
        }
    }


    
    delete A;
    delete R;

    cudaFree(d_D);
    cudaFree(d_G);
    cudaFree(d_Q);
    cudaFree(d_A);
    cudaFree(d_R);

    while(!PQ.empty()){
         auto x=PQ.top();
         V.push_back(x);
         //cout<<x.first<<" "<<get<3>(x.second)<<" "<<get<0>(x.second)<<" "<<get<1>(x.second)<<" "<<get<2>(x.second)<<"\n";
         PQ.pop();
    }

    ofstream fout;
    fout.open("output.txt");
    for(int i=V.size()-1;i>=0;i--){
        auto x=V[i];
        fout<<get<0>(x.second)<<" "<<get<1>(x.second)<<" "<<get<2>(x.second)<<"\n";
    }

    fout.close();

    auto end=std::chrono::high_resolution_clock::now();
    
    cout<<"Time taken for execution :"<<(1e-6*(std::chrono::duration_cast<std::chrono::nanoseconds>(end-beg)).count())<<" ms"<<"\n";
    

   
   

}