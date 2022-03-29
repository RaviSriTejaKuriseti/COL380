#include "psort.h"
#include <omp.h>


uint32_t minimum(uint32_t a,uint32_t b){
    if(a>=b){
        return b;
    }
    else{
        return a;
    }
}




void InsertionSort(uint32_t *arr,uint32_t left_bound,uint32_t right_bound){

    if(left_bound!=0){

         for(uint32_t i=left_bound+1;i<=right_bound;i++){

            uint32_t temp=arr[i];
            uint32_t j=i-1;
            while(j>=left_bound && arr[j]>temp){
                arr[j+1]=arr[j];
                j--;
            }
            arr[j+1]=temp;

        }

    }
    else{

        for(int i=1;i<=right_bound;i++){
            uint32_t temp=arr[i];
            int j=i-1;
            while(j>=0 && arr[j]>temp){
                arr[j+1]=arr[j];
                j--;
            }
            arr[j+1]=temp;

        }
    }

   


}


void Merge(uint32_t *arr,uint32_t start,uint32_t mid,uint32_t end){
            

        uint32_t l1=mid-start+1;
        uint32_t l2=end-mid;


        uint32_t *A=new uint32_t [l1];
        uint32_t *B=new uint32_t [l2];

        for(auto i=0;i<l1;i++){
            A[i]=arr[start+i];
        }
        for(auto I=0;I<l2;I++){
            B[I]=arr[mid+1+I];
        }

        uint32_t j=0;
        uint32_t k=0;
        uint32_t u=start;

        while(j<l1 && k<l2){
            if(A[j]<=B[k]){
                arr[u]=A[j];
                u++;
                j++;
            }
            else{
                arr[u]=B[k];
                u++;
                k++;
            }
        }




        while(j<l1){
            arr[u]=A[j];
            u++;
            j++;
        }

        while(k<l2){
            arr[u]=B[k];
            k++;
            u++;
        }


        delete A;
        delete B;

      


    
}


void HybridSort(uint32_t *arr,uint32_t length){  //Hybrid of Mergesort and InsertionSort.

    uint32_t RUN=0;
    if(length<=(1<<8)){
        RUN=16;
    }
    else if(length<=(1<<16)){
        RUN=32;
    }
    else{
        RUN=64;
    }

    for(uint32_t i=0;i<length;i+=RUN){
        InsertionSort(arr,i,minimum(i+RUN-1,length-1));
    }
    


    for(uint32_t j=RUN;j<length;j=2*j){

        for(uint32_t k=0;k<length;k+=2*j){

            uint32_t left=k;
            uint32_t mid=k+j-1;
            uint32_t right=minimum(2*j-1+k,length-1);
            if(mid<right){
                Merge(arr,left,mid,right);
            }

           
            

        }

    }

    
}





void ParallelSort(uint32_t *data, uint32_t n, int p)

{
    // Entry point to your sorting implementation.
    // Sorted array should be present at location pointed to by data.

    int thresh=(2*n)/p;
    uint32_t *left_ends = new uint32_t[p];  //[0,a1],[a1+1,a2],[a2+1,a3],.........[a(p-1)+1,a(p)]
    int rem=(n%p);
    int q=n/p;
    for(int i=0;i<p;i++){
        if(i<rem){            
            left_ends[i]=(q+1)*(i);
        }
        else{
            left_ends[i]=(q)*i+rem;

        }
        
    }




    uint32_t *pseudo_splitters= new uint32_t[p*p];
    for(int u=0;u<p;u++){
        for(int v=0;v<p;v++){
            pseudo_splitters[u*p+v]=data[left_ends[u]+v];
        }
    }

    delete left_ends;

        
    

    HybridSort(pseudo_splitters,p*p);

       
    

    uint32_t *Splitters= new uint32_t[p-1];
   


    for(int I=0;I<p-1;I++){
        Splitters[I]= pseudo_splitters[(I+1)*p];
    }

    

   
    delete pseudo_splitters;

    uint32_t *Sizes = new uint32_t[p];
    uint32_t *Start_Indices = new uint32_t[p];
    uint32_t *Ctr = new uint32_t[p];
    uint32_t *arr = new uint32_t[n];


    for(uint32_t I=0;I<p;I++){
        Sizes[I]=0;
    }


   
    for(uint32_t I=0;I<n;I++){
        if(data[I]<=Splitters[0]){
                Sizes[0]+=1;
            }
            else if(data[I]>Splitters[p-2]){
                Sizes[p-1]+=1;
            }
            else{
                for(int J=1;J<p-1;J++){
                    
                    if(data[I]>Splitters[J-1] && data[I]<=Splitters[J]){
                        Sizes[J]+=1;
                    }
                    else{
                        Sizes[J]+=0;
                    }
            
                }



            }
        
       
    }

   


    Start_Indices[0]=0;
    for(int I=0;I<p-1;I++){
        Start_Indices[I+1]=Start_Indices[I]+Sizes[I];
    }

   


   


    for(int I=0;I<p;I++){
       Ctr[I]=Start_Indices[I];
    }


    for(uint32_t I=0;I<n;I++){

        if(data[I]<=Splitters[0]){
            arr[Ctr[0]]=data[I];
            Ctr[0]++;
        }
        else if(data[I]>Splitters[p-2]){
            arr[Ctr[p-1]]=data[I];
            Ctr[p-1]++;
        }
        else{
            for(int J=1;J<p-1;J++){                
        
                if(data[I]>Splitters[J-1] && data[I]<=Splitters[J]){
                    arr[Ctr[J]]=data[I];
                    Ctr[J]++;
                }
            }
      
        }

    }

    delete Ctr;

   


      

 #pragma omp parallel
 
 
 
 {  
   
        
    for(int K=0;K<p;K++){

        #pragma omp task
        {

            uint32_t temp=Start_Indices[K];
            uint32_t temp1=Sizes[K];
            uint32_t *B=new uint32_t[temp1];
            
                

            for(uint32_t L=temp;L<temp+temp1;L++){
                B[L-temp]=arr[L];

            }

            
            if(Sizes[K]<=thresh){
            HybridSort(B,temp1);
            }
            else{
                ParallelSort(B,temp1,p);
                
            }


            for(uint32_t L=temp;L<temp+temp1;L++){
                data[L]=B[L-temp];

            }

            delete B;

        }

    }

    
 }


 


 


    delete arr;
    delete Splitters;
    delete Start_Indices;
    delete Sizes;
    


}
    

    

