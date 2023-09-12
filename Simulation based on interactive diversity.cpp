// Prisoners Dilemma game on a small-world graph constructed from a square lattice 
// Some players are blocked to give their strategy (other players cannot adopt their strategy)

// standard include
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
//#include <windows.h>
using namespace std;
// define priority classes
#define NORMAL_PRIORITY_CLASS       0x00000020
#define IDLE_PRIORITY_CLASS         0x00000040
#define HIGH_PRIORITY_CLASS         0x00000080
#define REALTIME_PRIORITY_CLASS     0x00000100

// define parameters
#define L         300      /* lattice size                   */
#define SIZE        (L*L)    /* number of sites                */
#define MC_STEPS    40000   /* run-time in MCS     */
#define Last_STEPS    10000
//#define r               /* temptation to defect */
#define K           0.1      /* temperature */
#define Q           0      /* Q portion of links are rewired */
#define NAMEOUT     "K4b075r5Q2"
#define RANDOMIZE   3145215
#define deta        0.3


int sumN,sumL,cooperator,defector,loner,cooperator0,defector0,loner0,cooperator1,defector1,loner1;
double b;
double u;   /* 边玩家的比例    */
double sigma=0.3;     /* 矩阵中loner的收益    */
double rho;  /* 以rho的概率学习其他邻居的策略，（1-rho）的概率学习直接相连的邻居的策略*/
double delta=0.000001;  /* 突变率                  */



typedef int       tomb1[SIZE];
typedef long int  tomb3[SIZE][4];
typedef double    tomb4[SIZE];


tomb1 player_s;           /* matrix ,containing players strategies */
tomb3 player_n;           /* matrix, containing players neighbours */
tomb1 player_Type;
tomb3 player_ls;
tomb1 player_ns;

void prodgraph(void);      /* creates host graph                    */
void initial(void);        /* initial state                         */
void game(void);
void tongji(void);

ofstream outfile1;
ofstream outfile2;


//以下是随机数产生模块，不用管它,直接用就行，用randf()可以直接产生0-1满足均匀分布的随机数，randi(x),产生0---x-1的随机整数
/*************************** RNG procedures ****************************************/
#define NN 624
#define MM 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long mt[NN]; /* the array for the state vector  */
static int mti=NN+1; /* mti==NN+1 means mt[NN] is not initialized */
void sgenrand(unsigned long seed)
{int i;
 for (i=0;i<NN;i++) {mt[i] = seed & 0xffff0000; seed = 69069 * seed + 1;
                     mt[i] |= (seed & 0xffff0000) >> 16; seed = 69069 * seed + 1;
  }
  mti = NN;
}
void lsgenrand(unsigned long seed_array[])
{ int i; for (i=0;i<NN;i++) mt[i] = seed_array[i]; mti=NN; }
double genrand() 
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    if (mti >= NN) 
    {
        int kk;
        if (mti == NN+1) sgenrand(4357); 
        for (kk=0;kk<NN-MM;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MM] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<NN-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(MM-NN)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[NN-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[NN-1] = mt[MM-1] ^ (y >> 1) ^ mag01[y & 0x1];
        mti = 0;
    }  
    y = mt[mti++]; y ^= TEMPERING_SHIFT_U(y); y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C; y ^= TEMPERING_SHIFT_L(y);
    return y;  
}

double randf(){ return ( (double)genrand() * 2.3283064370807974e-10 ); }
long randi(unsigned long LIM){ return((unsigned long)genrand() % LIM); }

/********************** END of RNG ************************************/


void initial(void)
{
	 int i,j,k;
     sumN=0;
	 sumL=0; 
	 cooperator=0;
	 defector=0;
	 loner=0;
	 cooperator0=0;
	 defector0=0;
	 loner0=0;
	 cooperator1=0;
	 defector1=0;
	 loner1=0;
	 
    for (i=0; i<SIZE; i++)
	{ 

		player_Type[i]=(int)randi(2);
		if(randf()<u)
		{
			player_Type[i]=1;  //边玩家 
		}
		else
		{
			player_Type[i]=0;  //点玩家 
		}
		if(player_Type[i]==0)
		{
		   sumN++;
		   player_ns[i]=(int)randi(3);
		   if(player_ns[i]==0) 
		   {
		   	 cooperator0=cooperator0+4;  //点策略合作者边数  
		   	 cooperator=cooperator+4;
		   }
		   else if(player_ns[i]==1)
		   {
		     defector0=defector0+4;
		     defector=defector+4;
		   }
		   else
		   {
		   	  loner0=loner0+4;
		   	  loner=loner+4;
		   }
		}
		else
		{
			sumL++;
		   for(j=0;j<4;j++)
		   {
		   	 player_ls[i][j]=(int)randi(3);
		   	 if(player_ls[i][j]==0)
		   	 { 	
				cooperator1++;	
				cooperator++;	
			 }
			 else if(player_ls[i][j]==1)
			 {
			    defector1++;	
			    defector++;	
			 }
			 else
			 {
			 	loner1++;	
				loner++;
			 }
		   }     
		}	
	}
}



// creates first a square grid graph and then rewires Q links 
void prodgraph(void)             
{
	int nneighbor, iu, ju, neighbor1, neighbor2;
	long int rewire, first, player1,player2,player3,MCe;
	double x;
	int i,j;
   

	// set up an initial square lattice
	for(i=0; i<L; i++)                     
	{
		for(j=0; j<L; j++)
		{ 
			// the first player
			player1 = L * j + i;             

			// and its four nearest neighbors
			iu = i + 1;         
			ju = j;     
			if (iu==L) iu = 0;
			player2 = L * ju + iu;  
			player_n[player1][0] = player2;

			iu = i;             
			ju = j + 1; 
			if (ju==L) ju = 0;
			player2 = L * ju + iu;  
			player_n[player1][1] = player2;

			iu = i - 1;         
			ju = j;     
			if (i==0) iu = L - 1;
			player2 = L * ju + iu;  
			player_n[player1][2] = player2;

			iu = i;             
			ju = j - 1; 
			if (j==0) ju = L - 1;
			player2 = L * ju + iu;  
			player_n[player1][3] = player2;
		}
	}

	// the rewiring starts - Q portion of joints is chosen randomly
	x = (double) (1.0 - Q);
	x = - log ( x );
	x = x*(2*SIZE);
	MCe = (long int) x;
	

	first = randi(SIZE);
	nneighbor = randi(4);
	player1 = player_n[first][nneighbor];
	neighbor1 = (nneighbor+2) % 4;

	for(rewire=0; rewire<MCe-1; rewire++)
	{
		do
		player2 = randi(SIZE);
		while(player2==first || player2==player1);
		do
		{
		neighbor2 = randi(4);
		player3 = player_n[player2][neighbor2];
		}
		while(player3==first || player3==player1);
		player_n[player1][neighbor1] = player2;
		player_n[player2][neighbor2] = player1;
		player1 = player3;
		if (player_n[player3][0]==player2)   neighbor1=0;
		if (player_n[player3][1]==player2)   neighbor1=1;
		if (player_n[player3][2]==player2)   neighbor1=2;
		if (player_n[player3][3]==player2)   neighbor1=3;
	 }	

	player_n[player1][neighbor1] = first;
	player_n[first][nneighbor] = player1;
  
	//cout<<player1<<'\t'<<player_n[player1][neighbor1]<<endl;
	
}

void game(void)
{
    int i,j,k,lo,h;
	int strat1,strat2;
	int player1,player2;
	double P1,P2,F1,F2;
    int player11,player22;
    int strat11,strat22;
    int suiji;                 //随机数
	double p,dP;

/*	ofstream outfile3("细分.txt",ios::out|ios::app);

		 if(!outfile3)
		{
        cout<<"不能打开此文件!";
        exit(1);
		}     */     
            

	    for (i=0; i<SIZE; i++)
		{
			player1 = (int) randi(SIZE);  
			   
		    P1=0;
			F1=0;
			P2=0;
			F2=0;
		//	strat1 = player_s[player1];
			
			if(player_Type[player1]==0)
			{
				strat1 = player_ns[player1];
				if(strat1==0)
				{
					for(j=0;j<4;j++)
					{
						player11=player_n[player1][j];
						if(player_Type[player11]==0)
						{
						   strat11=player_ns[player11];
						}
				        else
				        {
				           if(j<2)  lo=j+2;
				           else     lo=j-2;
						   strat11=player_ls[player11][lo];
						}
						
						if(strat11==0)        P1+=1; 
						else if(strat11==1)   P1+=0;	
						else                  P1+=sigma;
					}
				}
				else if(strat1==1)
				{
				    for(j=0;j<4;j++)
					{
						player11=player_n[player1][j];
						if(player_Type[player11]==0)
						{
						   strat11=player_ns[player11];
						}
				        else
				        {
				           if(j<2)  lo=j+2;
				           else     lo=j-2;
						   strat11=player_ls[player11][lo];
						}
						
						if(strat11==0)         P1+=b; 
						else if(strat11==1)    P1+=0;
						else	               P1+=sigma;
					}	
				}
				else
				{
					P1+=4*sigma;	 
				}
				F1=P1;    //点类型中心个体收益 
				
				suiji=(int)randi(4);   
				player2 = player_n[player1][suiji];  
                if (player_Type[player2]==0) 
                {
                	  strat2 = player_ns[player2];
				      if(strat2==0)
				      {
					      for(j=0;j<4;j++)
					      {
						     player22=player_n[player2][j];
						     if(player_Type[player22]==0)
						     {
						         strat22=player_ns[player22];
						     }
				             else
				             {
				                 if(j<2)  lo=j+2;
				                 else     lo=j-2; 
						         strat22=player_ls[player22][lo];
						     }
						     if(strat22==0)          P2+=1; 
						     else if(strat22==1)     P2+=0;
							 else 	                 P2+=sigma; 
					       }
				       }
				       else if(strat2==1)
				       {
				            for(j=0;j<4;j++)
					        {
						        player22=player_n[player2][j];
						        if(player_Type[player22]==0)
						        {
						           strat22=player_ns[player22];
						        }
				                else
				                {
				                    if(j<2)  lo=j+2;
				                    else     lo=j-2;
						            strat22=player_ls[player22][lo];
						        }
						
						        if(strat22==0)        P2+=b; 
						        else if(strat22==1)   P2+=0;
								else                  P2+=sigma;	
				        	}	
				        }
				        else
				        {
				        	P2+=4*sigma;	
						}
				}
			    else
			    {
			    	    for(j=0;j<4;j++)
			         	{
				  	        strat2 = player_ls[player2][j];
				  	        player22=player_n[player2][j];
				  	        if(player_Type[player22]==0)
				  	        {
				  	            strat22=player_ns[player22];	
					        }
					        else
					        {
						        if(j<2)  lo=j+2;
				                else     lo=j-2;
				  	            strat22=player_ls[player22][lo];
					        }
				  	
				  	        if(strat2==0)
				  	        {
				  	 	        if(strat22==0)          P2+=1; 
						        else if(strat22==1)     P2+=0;
						        else                    P2+=sigma;
					        }
					        else if(strat2==1)
					        {
					            if(strat22==0)          P2+=b; 
						        else if(strat22==1)     P2+=0;	
						        else                    P2+=sigma;
					        } 
							else
							{
								P2+=sigma;
							}   
				        }
				        if(suiji<2)  lo=suiji+2;
				        else         lo=suiji-2;
						strat2 = player_ls[player2][lo];
				}
			    F2=P2;       //邻居收益 
			    
			    
				if(randf()<delta)
				{
					player_ns[player1]=(int)randi(3);
				}
				else
				{
					if(strat1!=strat2)
			        {
					  dP=F1-F2;
			          p=1/(1+exp(dP/K));
			          if(randf()<p)
					  {
					  	player_ns[player1]=strat2;
					  } 
			        }              //策略更新
				}
			}
			else     //边类型收益计算 
			{
				for(j=0;j<4;j++)
				{
				  	strat1 = player_ls[player1][j];
				  	player11=player_n[player1][j];
				  	if(player_Type[player11]==0)
				  	{
				  	    strat11=player_ns[player11];	
					}
					else
					{
						if(j<2)  lo=j+2;
				        else     lo=j-2;
				  	    strat11=player_ls[player11][lo];
					}
				  	
				  	if(strat1==0)
				  	{
				  		if(strat11==0)           P1+=1; 
						else if(strat11==1)      P1+=0;
						else                     P1+=sigma;
					}
					else if(strat1==1)
					{
					    if(strat11==0)            P1+=b; 
						else if(strat11==1)       P1+=0;
						else                      P1+=sigma;	
					}
					else
					{
						P1+=sigma;	
					}
				}
				F1=P1;
				
				for(k=0;k<4;k++)
				{
					P2=0;
					strat1 = player_ls[player1][k];
					
					if(randf()<rho) 
					{
						h=(int)randi(3);
						if(k==0) 
						{
							h=h+1;
						}
						else if(k==1)
						{
							h=h+2;
							if(h>3)  h=h-4;
						}
						else if(k==2)
						{
							h=h+3;
							if(h>3)  h=h-4;
	                    }
	                    else
	                    {
	                    	h=h+0; 
						} 
					}
					else
					{
						h=k;
					}
					
					player2 = player_n[player1][h]; 
					
					if(player_Type[player2]==0)
					{
						strat2 = player_ns[player2];
				        if(strat2==0)
				        {
					        for(j=0;j<4;j++)
					        {
						         player22=player_n[player2][j];
						         if(player_Type[player22]==0)
						         {
						            strat22=player_ns[player22];
						         }
				                 else
				                 {
				                    if(j<2)  lo=j+2;
				                    else     lo=j-2; 
						            strat22=player_ls[player22][lo];
						         }
						         if(strat22==0)         P2+=1; 
						         else if(strat22==1)    P2+=0;
								 else	                P2+=sigma;
					        }
				        }
				        else if(strat2==1)
				        {
				            for(j=0;j<4;j++)
					        {
						        player22=player_n[player2][j];
						        if(player_Type[player22]==0)
						        {
						           strat22=player_ns[player22];
						        }
				                else
				                {
				                    if(j<2)  lo=j+2;
				                    else     lo=j-2;
						            strat22=player_ls[player22][lo];
						        }
						
						        if(strat22==0)           P2+=b; 
						        else if(strat22==1)      P2+=0;	
						        else                     P2+=sigma;	
				        	}	
				        }
				        else
				        {
				        	P2+=4*sigma;	
						}
					}
					else
					{
					        for(j=0;j<4;j++)
			         	    {
				  	            strat2 = player_ls[player2][j];
				  	            player22=player_n[player2][j];
				  	            if(player_Type[player22]==0)
				  	            {
				  	                strat22=player_ns[player22];	
					            }
					            else
					            {
						            if(j<2)  lo=j+2;
				                    else     lo=j-2;
				  	                strat22=player_ls[player22][lo];
					            }
				  	
				  	            if(strat2==0)
				  	            {
				  	 	            if(strat22==0)          P2+=1; 
						            else if(strat22==1)     P2+=0;
						            else                    P2+=sigma;
					            }
					            else if(strat2==1)
					            {
					                if(strat22==0)           P2+=b; 
						            else if(strat22==1)      P2+=0;	
						            else                     P2+=sigma;
					            }
								else
								{
									P2+=sigma;
								}    
				            }
				            if(h<2)      lo=h+2;
				            else         lo=h-2;
						    strat2 = player_ls[player2][lo];	
					} 
					F2=P2; 
					
					if(randf()<delta)
					{
						player_ls[player1][k]=(int)randi(3);
					}
					else
					{
						if(strat1!=strat2)
		                {
				           dP=F1-F2;
		                   p=1/(1+exp(dP/K));
		                   if(randf()<p)
				           {
				       	        player_ls[player1][k]=strat2;
				           }     
					    }
					}
				}    //一圈玩家2的收益 
			}


//*****************************************************随机选择邻居
 	 
 /*  static long mm = 1;
	           if(mm%800==0) {
		           cooperator=defector=loner=0;

		           for(int ii=0;ii<SIZE;ii++)
						{
	                        if(player_s[ii]==0)   cooperator++;
				            if(player_s[ii]==1) defector++;
                            if(player_s[ii]==2) loner++;
						}

			      double m,xxx,YYY,ZZZ;
		          m=(double)mm/SIZE;
			      xxx=(double)cooperator/SIZE;
				  YYY=(double)defector/SIZE;
				  ZZZ=(double)loner/SIZE;

	    //cout<<m<<'\t'<<xxx<<'\t'<<endl;
		          outfile3<<m<<'\t'<<ZZZ<<'\t'<<YYY<<'\t'<<xxx<<endl;
			   } 
	           mm += 1;
	//if(mm==(M*M+1)) mm=1;  */

           
		} //end of the SIZE  

}

           
void tongji(void)
{
     int i,j;

     sumN=0;
	 sumL=0; 
	 cooperator=0;
	 defector=0;
	 loner=0;
	 cooperator0=0;
	 defector0=0;
	 loner0=0;
	 cooperator1=0;
	 defector1=0;
	 loner1=0;
  

  for(i=0;i<SIZE;i++)
  {
        if(player_Type[i]==0)
		{
		   sumN++;
		   if(player_ns[i]==0) 
		   {
		   	 cooperator0=cooperator0+4;  //点策略合作者边数 
			 cooperator=cooperator+4; 
		   }
		   else if(player_ns[i]==1)
		   {
		     defector0=defector0+4;
		     defector=defector+4;
		   }
		   else
		   {
		   	 loner0=loner0+4;
		   	 loner=loner+4;
		   }
		}
		else
		{
			sumL++;
		   for(j=0;j<4;j++)
		   {
			   	 if(player_ls[i][j]==0)
			   	 { 	
					cooperator1++;
					cooperator++;	
				 }
				 else if(player_ls[i][j]==1)
				 {
				    defector1++;
					defector++;	
				 }
				 else
				 {
				 	loner1++;
				 	loner++;
				 }
		   }     
		}
  }

}
// the main program
int main()
{
	int steps;
    double t0,t1,aa,bb,cc,aa0,bb0,cc0,aa1,bb1,cc1,fN,fL,x0,x1,x2,xx0,xx1,xx2,xxx0,xxx1,xxx2,TN,TL,XC,XD,XL,XC0,XD0,XL0,XC1,XD1,XL1;
    

	outfile1.open("frequency.txt");
	outfile2.open("average.txt");

    if(!outfile1)
	{
	 cout<<"can not open";
	 abort();
	}
   
	if(!outfile2)
	{
	 cout<<"can not open";
	 abort();
	}

	// initialize the random number generation
	sgenrand(RANDOMIZE);

	prodgraph();

	



	// begins the mc steps
for(rho=0.0;rho<=1.001;rho=rho+0.025)
{
 for(u=0.0;u<=1.01;u=u+0.025)
 {
  for(b=1.0;b<=2.001;b=b+0.025)
  {
    t0=0;
    t1=0;
	
	aa=0;
	bb=0;
	cc=0;
	
	aa0=0;
	bb0=0;
	cc0=0;
	
	aa1=0;
	bb1=0;
	cc1=0;
	
	initial();
	for (steps=0; steps<MC_STEPS; steps++)
	{
	 game();
	 tongji();

	 fN=(double)sumN/SIZE;
     fL=(double)sumL/SIZE;
     
     x0=(double)cooperator/(4*SIZE);
     x1=(double)defector/(4*SIZE);
     x2=(double)loner/(4*SIZE);
     
	 xx0=(double)cooperator0/(4*SIZE);
     xx1=(double)defector0/(4*SIZE);
     xx2=(double)loner0/(4*SIZE);
     
     xxx0=(double)cooperator1/(4*SIZE);
     xxx1=(double)defector1/(4*SIZE);
     xxx2=(double)loner1/(4*SIZE);
     
	 outfile1<<steps<<'\t'<<delta<<'\t'<<rho<<'\t'<<u<<'\t'<<b<<'\t'<<fN<<'\t'<<fL<<'\t'<<x0<<'\t'<<x1<<'\t'<<x2<<'\t'<<xx0<<'\t'<<xx1<<'\t'<<xx2<<'\t'<<xxx0<<'\t'<<xxx1<<'\t'<<xxx2<<endl;

     
           if(steps>MC_STEPS-Last_STEPS-1)
		   {
		     t0+=fN;
		     t1+=fL;
		     
		     aa+=x0;
		     bb+=x1;
		     cc+=x2;
		     
		     aa0+=xx0;
		     bb0+=xx1;
		     cc0+=xx2;
		     
		     aa1+=xxx0;
		     bb1+=xxx1;
		     cc1+=xxx2;
		     
		     TN=(double)t0/Last_STEPS;
		     TL=(double)t1/Last_STEPS;
		     
		     XC=(double)aa/Last_STEPS;
		     XD=(double)bb/Last_STEPS;
		     XL=(double)cc/Last_STEPS;
		     
		     XC0=(double)aa0/Last_STEPS;
		     XD0=(double)bb0/Last_STEPS;
		     XL0=(double)cc0/Last_STEPS;
		     
		     XC1=(double)aa1/Last_STEPS;
		     XD1=(double)bb1/Last_STEPS;
		     XL1=(double)cc1/Last_STEPS;
		   }
	}//end of the MC steps
    outfile2<<delta<<'\t'<<rho<<'\t'<<u<<'\t'<<b<<'\t'<<TN<<'\t'<<TL<<'\t'<<XC<<'\t'<<XD<<'\t'<<XL<<'\t'<<XC0<<'\t'<<XD0<<'\t'<<XL0<<'\t'<<XC1<<'\t'<<XD1<<'\t'<<XL1<<endl;
	cout<<delta<<'\t'<<rho<<'\t'<<u<<'\t'<<b<<'\t'<<TN<<'\t'<<TL<<'\t'<<XC<<'\t'<<XD<<'\t'<<XL<<'\t'<<XC0<<'\t'<<XD0<<'\t'<<XL0<<'\t'<<XC1<<'\t'<<XD1<<'\t'<<XL1<<endl;
  
   }
  }
}
	outfile1.close();
	outfile2.close();
	
	return 0;
}	
