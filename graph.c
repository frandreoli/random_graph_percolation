#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define N_0 100
#define SEED 1651265127
#define F_FAILURE -1
#define F_SUCCESS 0
#define N_GRAPH 20
#define DELETE
#define WIN_OS 0


//DICHIARAZIONI STRUTTURE
struct connection{unsigned long long int n_bridges; int *bridges;};
struct graph {struct graph *previous; unsigned long long int size; struct connection neigh;} *nodes;
struct solution {double mean; double std_dev; double value[N_GRAPH];} read_smax, read_ssqmean;


//VARIABILI GLOBALI
const unsigned long long int N=N_0;
unsigned long long int M=0;
unsigned long long int size_max;
unsigned long long int size_square_mean;
double c=0.;


//FUNZIONI
void error(char *);
void error(char *message){
printf("\n\n%s.\n\n",message);
exit(EXIT_FAILURE);}

unsigned long long int uli_random(unsigned long long int, unsigned long long int);
unsigned long long int uli_random(unsigned long long int a, unsigned long long int b){
return a+((unsigned long long int)rand())%(b-a);}

void initialization(void);
void initialization(void){
int i;
M=0;
size_max=1;
size_square_mean=N; 
for(i=0; i<N; i++){
nodes[i].previous=&(nodes[i]);
nodes[i].size=1;
nodes[i].neigh.n_bridges=0;}}

int check_bridge(unsigned long long int, unsigned long long int);
int check_bridge(unsigned long long int node1, unsigned long long int node2){
int i;
for(i=0;i<nodes[node1].neigh.n_bridges;i++){
if(nodes[node1].neigh.bridges[i]==node2) return F_FAILURE;}
return F_SUCCESS;}

void add_bridge(unsigned long long int, unsigned long long int);
void add_bridge(unsigned long long int node1, unsigned long long int node2){
int n_b1=++(nodes[node1].neigh.n_bridges);
int n_b2=++(nodes[node2].neigh.n_bridges);
nodes[node1].neigh.bridges=realloc(nodes[node1].neigh.bridges, n_b1*sizeof(int));
nodes[node2].neigh.bridges=realloc(nodes[node2].neigh.bridges, n_b2*sizeof(int));
nodes[node1].neigh.bridges[n_b1-1]=node2;
nodes[node2].neigh.bridges[n_b2-1]=node1;}

struct graph *head(struct graph *);
struct graph *head(struct graph *node){
while(node->previous!=node) {node = node->previous;}
return node;}

void cluster_develop(unsigned long long int, unsigned long long int);
void cluster_develop(unsigned long long int node1, unsigned long long int node2){
unsigned long long int new_size;
struct graph *head1=head(&(nodes[node1]));
struct graph *head2=head(&(nodes[node2]));
M++;
if(head1!=head2){
new_size = head1->size + head2->size;
if(new_size>size_max) size_max=new_size;
size_square_mean+=2*(unsigned long long int)head1->size*head2->size;     
if(head1->size>=head2->size) {head2->previous=head1; head1->size=new_size;}
else {head1->previous=head2; head2->size=new_size;}
}}

void step(void);
void step(void){
unsigned long long int node1=uli_random(0,N), node2;
do{node2=uli_random(0,N);} while(node1==node2 || check_bridge(node1,node2)==F_FAILURE);
add_bridge(node1,node2);
cluster_develop(node1,node2);
c=2.*((double)M/N);}

double std_dev_calc(struct solution);
double std_dev_calc(struct solution sol){
double std_dev=0;
int i;
for(i=0;i<N_GRAPH;std_dev+=(sol.value[i]-sol.mean)*(sol.value[i++]-sol.mean));
std_dev=sqrtl(std_dev/(double)(N_GRAPH-1.));
return std_dev/sqrtl((double)(N_GRAPH)); //dev std sulla media!!
}



//MAIN
void main(){
int i=0, j=0, end_file;
FILE **fp_array, *fp_single;
char f_string[20], temp;
#ifdef SEED
srand(SEED);
#endif
#ifndef SEED
srand(time(0));
#endif
char rm_string[4];
if(WIN_OS)
    {sprintf(rm_string, "del");}
else
    {sprintf(rm_string, "rm");}

//Alloco memoria
nodes=(struct graph *)calloc(N, sizeof(struct graph));
fp_array=(FILE **)calloc(N_GRAPH, sizeof(FILE *));

//Calcolo traiettorie
for(i=0;i<N_GRAPH;i++){
initialization();
sprintf(f_string,"Data//graph%d.dat",i);
if((fp_single=fopen(f_string,"w+"))==NULL) error("ERROR: CANNOT OPEN DATA FILE");
fprintf(fp_single,"#Smax       S^2mean       c\n");
for(c=0;c<2;){step();
fprintf(fp_single,"%g        %g        %g\n", (double)size_max/N, (double)(size_square_mean-(unsigned long long int)size_max*size_max)/N, c);}
fclose(fp_single);}

printf("\n\nHERE I AMMMMMMMMMMMMMMMM");

//Imposto puntatori e buffer di lettura
for(i=0;i<N_GRAPH;i++){
sprintf(f_string,"Data//graph%d.dat",i);
if((fp_array[i]=fopen(f_string,"r"))==NULL) 
{printf("\n\nFILE ERROR: %s",f_string); error("ERROR: CANNOT OPEN-READ DATA FILE");}
}
sprintf(f_string,"g_mean_%Ld.dat",N);
if((fp_single=fopen(f_string, "w+"))==NULL) error("ERROR: CANNOT OPEN-WRITE DATA FILE");
fprintf(fp_single,"#Smax_mean          <S^2>_mean          c           dev_Smax_mean       dev_<S^2>_mean\n");
for(i=0;i<N_GRAPH;i++){
do{fscanf(fp_array[i],"%c",&temp);}while(temp!='\n');}

//Leggo e calcolo la media, e la dev standard
read_smax.mean=0.; read_ssqmean.mean=0.; read_smax.std_dev=0.; read_ssqmean.std_dev=0.;
for(end_file=0; end_file!=EOF; read_smax.mean=0., read_ssqmean.mean=0., read_smax.std_dev=0., read_ssqmean.std_dev=0.){
/*scorro i file dati*/
end_file=fscanf(fp_array[0], "%Lf   %Lf    %Lf", &read_smax.value[0], &read_ssqmean.value[0], &c);
read_smax.mean+=read_smax.value[0];
read_ssqmean.mean+=read_ssqmean.value[0];
for(i=1;i<N_GRAPH;i++){
fscanf(fp_array[i], "%Lf   %Lf    %*f", &read_smax.value[i], &read_ssqmean.value[i]);
read_smax.mean+=read_smax.value[i];
read_ssqmean.mean+=read_ssqmean.value[i];}
/*calcolo media*/read_smax.mean/=(double)N_GRAPH;
read_ssqmean.mean/=(double)N_GRAPH;
/*calcolo dev std della media*/read_smax.std_dev=std_dev_calc(read_smax);
read_ssqmean.std_dev=std_dev_calc(read_ssqmean);
if(end_file!=EOF)
fprintf(fp_single, "%Lg  %Lg  %Lg  %Lg  %Lg\n", read_smax.mean, read_ssqmean.mean, c, read_smax.std_dev, read_ssqmean.std_dev);}

//Cancello files
for(i=0;i<N_GRAPH;i++){
fclose(fp_array[i]);
#ifdef DELETE
sprintf(f_string,"%s Data//graph%d.dat",rm_string, i);
if(system(f_string)==-1) error("ERROR: CANNOT DELETE FILE");
#endif
}
fclose(fp_single);

}
