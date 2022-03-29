#include "classify.h"
#include <omp.h>

Data classify(Data &D, const Ranges &R, unsigned int numt)
{ // Classify each item in D into intervals (given by R). Finally, produce in D2 data sorted by interval
   assert(numt < MAXTHREADS);
   Counter counts[R.num()]; // I need on counter per interval. Each counter can keep pre-thread subcount.
   int N=D.ndata/numt;
   int rem=D.ndata%numt;
   int M=N+1;

   #pragma omp parallel num_threads(numt)
   {
      int tid = omp_get_thread_num(); // I am thread number tid
      // std::cout<<"Thread_num"<<" "<<tid<<"\n";

      if(tid<rem){


        


         for(int i=tid*M; i<(tid+1)*M; i++) { // Threads together share-loop through all of Data
            int v = D.data[i].value = R.range(D.data[i].key);// For each data, find the interval of data's key,
                          // and store the interval id in value. D is changed.
            counts[v].increase(tid); // Found one key in interval v
         }

      }

      else{




         for(int i=tid*N+rem; i<(tid+1)*N+rem; i++) { // Threads together share-loop through all of Data
            int v = D.data[i].value = R.range(D.data[i].key);// For each data, find the interval of data's key,
                          // and store the interval id in value. D is changed.
            counts[v].increase(tid); // Found one key in interval v
         }

      }

      


      
   }

   // Accumulate all sub-counts (in each interval;'s counter) into rangecount
   unsigned int *rangecount = new unsigned int[R.num()];
   for(int r=0; r<R.num(); r++) { // For all intervals
      rangecount[r] = 0;
      for(int t=0; t<numt; t++) // For all threads
         rangecount[r] += counts[r].get(t);
   //std::cout << rangecount[r] << " elements in Range " << r << "\n"; // Debugging statement
   }

   // Compute prefx sum on rangecount.
   for(int i=1; i<R.num(); i++) {
      rangecount[i] += rangecount[i-1];
   }

   // Now rangecount[i] has the number of elements in intervals before the ith interval.

   Data D2 = Data(D.ndata); // Make a copy

   int N1=R.num()/numt;
   int rem1=R.num()%numt;
   int M1=N1+1;
   
   #pragma omp parallel num_threads(numt)
   {
      int tid = omp_get_thread_num();
      //std::cout<<"Thread_num"<<" "<<tid<<"\n";

      if(tid<rem1){

         #pragma omp for

         for(int r=tid*M1; r<(tid+1)*M1; r++) { // Thread together share-loop through the intervals 
            int rcount = 0;
            for(int d=0; d<D.ndata; d++) // For each interval, thread loops through all of data and  
                if(D.data[d].value == r) // If the data item is in this interval 
                    D2.data[rangecount[r-1]+rcount++] = D.data[d]; // Copy it to the appropriate place in D2.
         }


      }

      else{

         #pragma omp for

          for(int r=tid*N1+rem1; r<(tid+1)*N1+rem1; r++) { // Thread together share-loop through the intervals 
            int rcount = 0;
            for(int d=0; d<D.ndata; d++) // For each interval, thread loops through all of data and  
                if(D.data[d].value == r) // If the data item is in this interval 
                    D2.data[rangecount[r-1]+rcount++] = D.data[d]; // Copy it to the appropriate place in D2.
         }


      }
      
         // for(int r=tid; r<R.num(); r+=numt) { // Thread together share-loop through the intervals 
         //    int rcount = 0;
         //    for(int d=0; d<D.ndata; d++) // For each interval, thread loops through all of data and  
         //        if(D.data[d].value == r) // If the data item is in this interval 
         //            D2.data[rangecount[r-1]+rcount++] = D.data[d]; // Copy it to the appropriate place in D2.
         // }

      
      
   }

   return D2;
}
