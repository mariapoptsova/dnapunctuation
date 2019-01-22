#include <err.h>		
#include <stdio.h>
#include <stdlib.h>
#include <string> 
#include <iostream>
#include <fstream>
#include <ctype.h>
#include <vector>
#include <sstream>
using namespace std;

FILE *fw;
int **M;
int minStem, maxStem, minLoop, maxLoop, maxGap;
int ROWS, COLUMNS, SLEN;
char *s1, *s2, *rez1, *rez2, *rez2Inverse;
int stemSize;
int exonInter;
string palStrAlnPrev;


/*
int pttSIZE;
int *startList, *endList, *startEx, *endEx, *exonCountList;
vector<string> giList; 
vector<string> exonStartList;
vector<string> exonEndList;
**/

      int MIN(int X, int Y){
          if(X<=Y)return X; else return Y;
       }

      int MAX(int X, int Y){
          if(X>=Y)return X; else return Y;
       }

      int Compl(char a, char b){
        if( ((a == 'T')&&(b == 'A')) || ((b == 'T')&&(a == 'A'))||
            ((a == 'C')&&(b == 'G')) || ((a == 'G')&&(b == 'C')) )
          return 1; else return 0; 
      }

// result is in rez2Inverse
      void inverse(char *s, int len){
         for(int i=0; i<len;i++){
              rez2Inverse[len-1-i]=s[i];
         }
      }

      char* shift(char *s, int len){
        for(int i=0; i<len-1;i++){
              s[i]=s[i+1];
         }
       return s;
      } 

      int printM()
      {
        int i,j;
        for(i=0;i<ROWS;i++){
         for(j=0;j<COLUMNS;j++){
           cout << M[i][j] << "\t";
         }
         cout << endl;
        } 
      }

/*
      int exonIntersect(int palStart, int palEnd, int index, int exonNu){
           int Intersect=0;
           stringstream ssStart(exonStartList[index]);
           stringstream ssEnd(exonEndList[index]);
           string field;
           startEx= new int[exonNu];
           endEx= new int[exonNu];
           //startExons
           int n=0;
            while (getline( ssStart, field, ',' ))
            {
               stringstream fs(field);
               int f;  //
               fs >> f;
               startEx[n]=f;
               n++;
            }
           //endExons
            n=0;
            while (getline( ssEnd, field, ',' ))
            {
               stringstream fs(field);
               int f;  //
               fs >> f;
               endEx[n]=f;
               n++;
            }
           //----------------------

           for(int i=0;i<exonNu;i++){
              if(((palStart>=startEx[i]&&palStart<=endEx[i])) ||
               ((palEnd>=startEx[i]&&palEnd<=endEx[i])))
                   Intersect=1;
           }
           delete startEx;
           delete endEx;
        return Intersect;
      }
 **/
 /*
      int geneIntersect(int palStart, int palEnd)
      { 
        exonInter=0;
        for(int i=0; i< pttSIZE; i++){
          if(((palStart>=startList[i]&&palStart<=endList[i])) || 
             ((palEnd>=startList[i]&&palEnd<=endList[i])))
           {
             /*
               if(exonCountList[i]>1){
                 exonInter=exonIntersect(palStart,palEnd,i, exonCountList[i]);
               }
             *//*
                 return i;
           }
          if(palEnd<startList[i]) break;         
        }
      return -1;
      }
 **/    

        void align(char *seq, int begin, int loop, int max_gap, int min_match, int max_match)
        {
            
//          string dna(seq);
//cout << dna << "This is dna\n";
           int i,j; 
            
           for (i = begin; i <= begin + max_match; i++ )
            {
                s1[max_match - i + begin] = seq[i];
            }
           
            for (i = 0; i < max_match + 1; i++ )
            {
                s2[i] = seq[begin + max_match + loop + 1+i];
            }
      

//string printS1(s1,SLEN);
//string printS2(s2,SLEN);
//cout << "s1=" << printS1 << "||ENDs1\n";
//cout << "s2=" << printS2 << "||ENDs2\n" ;
             
            int m = 0, w = 1;
            int low = 0, high = 0;

            for (i = 0; i <= SLEN; i++)
            {
                M[0][i] = i;
                M[i][0] = i;
            }

            for (i = 1; i <= SLEN; i++)
            {
                for (int j = 1; j <= SLEN; j++)
                    M[i][j] = 1000;
            }

            for (i = 1; i <= SLEN; i++)
            {
                low = 1;
                high = SLEN;

                for (j = low; j <= high; j++)
                {
                    if (Compl(s1[i - 1],s2[j - 1])) m = 0;
                    else m = 1;
                    M[i][j] = MIN((M[i - 1][j - 1] + m), MIN((M[i - 1][j] + w),(M[i][j - 1] + w)) );
                }
            }
//printM();
            int maxi = 0, maxj = 0, curi = 0, curj = 0;

            for (i = 1; i <= SLEN; i++)
            {
                for (j = 1; j <= SLEN; j++)
                {
                    if (M[i][j] <= max_gap && MIN(i, j) > MIN(maxi, maxj))
                    {
                        maxi = i;
                        maxj = j;
                    }
                }
            }

            if (MIN(maxi, maxj) >= min_match)
            {

                curi = maxi;
                curj = maxj;


                int k = 0;
                for (k = 0; curi > 0 && curj > 0; k++)
                {
                    if (M[curi][curj] == M[curi - 1][curj] + w)
                    {
                        rez1[k] = s1[curi - 1];
                        rez2[k] = '_';

                        curi = curi - 1;
                    }
                    else
                        if (M[curi][curj] == M[curi][curj - 1] + w)
                        {
                            rez1[k] = '_';
                            rez2[k] = s2[curj - 1];

                            curj = curj - 1;
                        }
                        else
                        {
                            rez1[k] = s1[curi - 1];
                            rez2[k] = s2[curj - 1];

                            curi = curi - 1;
                            curj = curj - 1;
                        }
                }

                for (; curi > 0; k++)
                {
                    rez1[k] = s1[curi - 1];
                    rez2[k] = '_';
                    curi = curi - 1;
                }

                for (; curj > 0; k++)
                {
                    rez1[k] = '_';
                    rez2[k] = s2[curj - 1];
                    curj = curj - 1;
                }
              
                stemSize=k;
                int k1 = 0, k2 = 0;
                for(i=0;i<k;i++)
                {
                  if(rez1[i] != '_') k1++;
                  if(rez2[i] != '_') k2++;
                  if(!Compl(rez1[i],rez2[i])){
                    rez1[i]=tolower(rez1[i]);
                    rez2[i]=tolower(rez2[i]);
                  }
                }

                int start, end;
                start=begin + max_match - k1+1;
                end=begin + max_match + loop + k2;
              // correction for ends
              // loop-end
              int loopActual=loop;
              if(islower(rez1[k-1]) || islower(rez2[k-1])){
                 loopActual=loopActual+2;
                 stemSize--;
                 k1--;
                 k2--;
               }  
              // stem-end
              if(islower(rez1[0]) || islower(rez2[0])){
                  rez1=shift(rez1,stemSize);
                  rez2=shift(rez2,stemSize);
                  k1--;
                  k2--;
                  stemSize--;
                  start++;
                  end--;
              }

            if((stemSize>=min_match)&&(stemSize<=max_match)&&(loopActual>=minLoop)&&(loopActual<=maxLoop))
            {
              // results are in rez2Inverse (errors withh memory in the previous versions)    
                inverse(rez2,stemSize);
                string rez1Str(rez1,stemSize); 
                string rez2Str(rez2Inverse,stemSize); 
             
                char *loopSeq = new char[loopActual];
                for(i=0;i<loopActual;i++){
                 loopSeq[i]=seq[start+k1+i];

                }
                string loopStr(loopSeq,loopActual);
                string palStrAln = rez1Str; 
                palStrAln+='-';
                if(loopActual!=0){
                    palStrAln+=loopStr;
                    palStrAln+='-';
                    palStrAln+=rez2Str;
                }else{
                    palStrAln+=rez2Str;
                }
  
                int palSize=k1+loopActual+k2;
                char *palSeq = new char[palSize];
                for(i=0;i<palSize;i++)
                { 
                   palSeq[i]=seq[start+i];
                }
                string palStr(palSeq,palSize);
                
                if(palStrAln != palStrAlnPrev)
                {
                    fprintf(fw,"%d\t%d\t%d\t%d\t",start,end,stemSize,loopActual);
                    fprintf(fw,"%s\t%s\t%s\t%s\t%s\t",rez1Str.c_str(),rez2Str.c_str(),loopStr.c_str(),palStrAln.c_str(),palStr.c_str());
					
					
				/* **/ fprintf(fw,"\n"); /* **/
				
				
				/*
                    int index;
                    index=geneIntersect(start,end);
                    if(index==-1){
                     fprintf(fw,"0\t");
                    }else{
                     fprintf(fw,"%s\t",giList[index].c_str());
                    } 

                   fprintf(fw,"%d\n",exonInter);
                   /*
                   cout << start << "\t" << end << "\t"<< stemSize<<"\t"<<loopActual<< "\t";
                   cout << rez1Str << "\t" << rez2Str << "\t" << loopStr << "\t" << palStrAln << "\t";
                   cout << palStr  << endl;
                  **/
                   palStrAlnPrev=palStrAln;
                    
                }
                delete palSeq; 
                delete loopSeq; 
             }               
             
            }
         }

int main(int argc, char* argv[])
{

	
	//size_t len, pttLen;
	string line, pttLine;
        int k=0;
        int i,size; 
        char *seq;

   if (argc < 7) {
          fprintf(stderr, "Usage: %s <genome.fna> <min_stem> <max_stem> <min_loop> <max_loop> <gap>\n", argv[0]);
          return 1;
        }


        minStem = atoi(argv[2]);
        maxStem = atoi(argv[3]);
        minLoop = atoi(argv[4]);
        maxLoop = atoi(argv[5]);
        maxGap = atoi(argv[6]);


        ifstream f( argv[1] ) ;
        if(!f.is_open())
               /**err(1, argv[1]);*/
				return 20;
        if (f) 
        {
 
          f.seekg (0, ios::end);
          size = f.tellg();
          f.seekg (0, ios::beg);
          seq = new char[size];
         int len;
	 while(getline(f, line)) 
         {
           if(line[0]!='>')
           {
              len=line.length();
              int i;
              for(i=0;i<len;i++){
               if((line[i]!='\n')&&(isalpha(line[i]))){
                 seq[k]=toupper(line[i]);
                 k++;
               }
              }  
            }
	  } 
        }
        f.close();
        int genomeSize=k;
 //fwrite(seq, genomeSize, 1, stdout);
 //printf("genomeSize=%d \n",genomeSize);
                     
     
/*
       //-------- ptt --------------------------
         ifstream fptt( argv[ 2 ] ) ;
         int LineNu=0;
     if(!fptt.is_open())
         err(1, argv[2]);
     if(fptt)
     {
         while(getline(fptt, pttLine))
         {
             LineNu++;
         }
        pttSIZE=LineNu;
 
//   cout << "pttSIZE: " << pttSIZE << endl;
        startList= new int[pttSIZE];
        endList = new int[pttSIZE];
        exonCountList = new int[pttSIZE];
        int start,end,cdsStart,cdsEnd,exonCount;
        char strand;
        char gi[12];
        char chr[6];
        char exonStarts[60000];  
        char exonEnds[60000];  
        int count=0;
        fptt.clear();
        fptt.seekg (0, ios::beg);
        while(getline(fptt, pttLine))
        {
              sscanf(pttLine.c_str(), "%s\t%s\t%c\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t", gi, chr,&strand, &start,&end,&cdsStart,&cdsEnd,&exonCount,exonStarts,exonEnds);
              startList[count]=start;
              endList[count]=end;
              exonCountList[count]=exonCount;
              exonStartList.push_back(exonStarts);
              exonEndList.push_back(exonEnds); 
              giList.push_back(gi);
              count++; 
        }
      }
       fptt.close();
      //-------- end of ptt ------------------
           
//cout  << "end of ptt" << endl; 
**/

     //----- create file --------------------
        char nf[250];
        char t[3]; 
        sscanf(argv[1],"%s.fna",nf);
        string nfStr(nf);
        nfStr+=".S";
        sprintf(t,"%d",minStem);  
        nfStr+=t;
        nfStr+="-";
        sprintf(t,"%d",maxStem);  
        nfStr+=t;
        nfStr+="_L";
        sprintf(t,"%d",minLoop);  
        nfStr+=t;
        nfStr+="-";
        sprintf(t,"%d",maxLoop);  
        nfStr+=t;
        nfStr+="_M";
        sprintf(t,"%d",maxGap);  
        nfStr+=t;
        nfStr+=".pal";
        
        fw = fopen(nfStr.c_str(), "w");
        if (fw == NULL)
                err(1, "nfStr.c_str()");


     //------------------------------------- 

            int l;
            k = genomeSize - maxStem - maxLoop - 1;
            M=0;
            SLEN=maxStem+1; 
            ROWS=SLEN+1;
            COLUMNS=ROWS;
             

            M= new int *[ROWS];
             for(i = 0 ; i < ROWS ; i++ )
               M[i] = new int[COLUMNS];

           s2 = new char[SLEN];            
           s1 = new char[SLEN];            

           rez1 = new char[SLEN + maxGap];
           rez2 = new char[SLEN + maxGap];
           rez2Inverse = new char[SLEN + maxGap];
 
           palStrAlnPrev=' ';

            for (i = maxStem; i < k; i++)
            {
              
                for (l = minLoop; l <= maxLoop; l++)
                {
                    align(seq, i - maxStem, l, maxGap, minStem, maxStem);
                }
            }

   delete seq;
   for(i = 0 ; i < ROWS ; i++ )
     delete [] M[i];
   delete [] M;
   delete [] s1;
   delete [] s2;
   delete [] rez1;
   delete [] rez2;
   delete [] rez2Inverse;
 /*
   delete [] startList;
   delete [] endList;
   delete [] exonCountList;
   **/
	return 0;
}
