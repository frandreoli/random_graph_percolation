#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define INT_TYPE unsigned /*long long*/ int
//#define N_NODES 100000
#define SEED 1651265127
#define F_FAILURE -1
#define F_SUCCESS 0
#define N_TRAJECTORIES 200
//Uncomment or comment the following line if you want to respectively delete or keep the trajectory files
#define DELETE
//Uncomment the following line if you work in Windows
#define WIN_OS


//Structure declaration
struct connection{INT_TYPE n_bridges; INT_TYPE *bridges;};
struct graph {struct graph *previous; INT_TYPE size; struct connection neigh;} *nodes;
struct solution {double mean; double std_dev; double value[N_TRAJECTORIES];} read_smax, read_ssqmean;


//Functions
void error(char *);
void error(char *message){
    printf("\n\n%s.\n\n",message);
    exit(EXIT_FAILURE);
}

INT_TYPE uli_random(INT_TYPE, INT_TYPE);
INT_TYPE uli_random(INT_TYPE a, INT_TYPE b){
    int i = 0;
    int resolution = 100;
    unsigned int final_drand=0;
    for(i=0;i<resolution;i++){
        final_drand+=rand();
    }
    return (INT_TYPE)(a+(final_drand)%(b-a));
}

void initialization(INT_TYPE);
void initialization(INT_TYPE N){
    int i;
    for(i=0; i<N; i++){
        nodes[i].previous=&(nodes[i]);
        nodes[i].size=1;
        nodes[i].neigh.n_bridges=0;
    }
}

int check_bridge(INT_TYPE, INT_TYPE);
int check_bridge(INT_TYPE node1, INT_TYPE node2){
    int i;
    for(i=0;i<nodes[node1].neigh.n_bridges;i++){
        if(nodes[node1].neigh.bridges[i]==node2) return F_FAILURE;
    }
    return F_SUCCESS;
}

void add_bridge(INT_TYPE, INT_TYPE);
void add_bridge(INT_TYPE node1, INT_TYPE node2){
    //Pre-incrementing ++() the number of neighbouring bridges 
    //and then assigning the value to the temporary variables n_b1 and n_b2
    INT_TYPE n_b1=++(nodes[node1].neigh.n_bridges);
    INT_TYPE n_b2=++(nodes[node2].neigh.n_bridges);
    //Adding memory to account for the new bridges
    nodes[node1].neigh.bridges=realloc(nodes[node1].neigh.bridges, n_b1*sizeof(INT_TYPE));
    nodes[node2].neigh.bridges=realloc(nodes[node2].neigh.bridges, n_b2*sizeof(INT_TYPE));
    //Adding the new bridges
    nodes[node1].neigh.bridges[n_b1-1]=node2;
    nodes[node2].neigh.bridges[n_b2-1]=node1;
}

struct graph *head(struct graph *);
struct graph *head(struct graph *node){
while(node->previous!=node) {
    node = node->previous;
    }
return node;
}

void cluster_develop(INT_TYPE, INT_TYPE, INT_TYPE, INT_TYPE *);
void cluster_develop(INT_TYPE node1, INT_TYPE node2, INT_TYPE M, INT_TYPE *size_array){
    INT_TYPE new_size;
    struct graph *head1=head(&(nodes[node1]));
    struct graph *head2=head(&(nodes[node2]));
    if(head1!=head2){
        new_size = head1->size + head2->size;
        if(new_size>size_array[0]) size_array[0]=new_size;
        size_array[1]+=2*(INT_TYPE)head1->size*head2->size;     
        if(head1->size>=head2->size) {
            head2->previous=head1; 
            head1->size=new_size;
            }
        else {
            head1->previous=head2; 
            head2->size=new_size;
            }
    }
}

void step(INT_TYPE, INT_TYPE, INT_TYPE *);
void step(INT_TYPE N, INT_TYPE M, INT_TYPE *size_array){
    //Stepp to add a new bridge between two random nodes that were not already connected
    INT_TYPE node1=uli_random(0,N), node2;
    do{
        node2=uli_random(0,N);
        } while(node1==node2 || check_bridge(node1,node2)==F_FAILURE);
    //Adding the bridge
    add_bridge(node1,node2);
    //Updating the cluster to account for the new bridge
    cluster_develop(node1,node2, M,size_array);
}

double std_dev_calc(struct solution);
double std_dev_calc(struct solution sol){
    double std_dev=0;
    int i;
    for(i=0;i<N_TRAJECTORIES;std_dev+=(sol.value[i]-sol.mean)*(sol.value[i++]-sol.mean));
    std_dev=sqrtl(std_dev/(double)(N_TRAJECTORIES-1.));
    return std_dev/sqrtl((double)(N_TRAJECTORIES)); 
}

INT_TYPE input_read(char *);
INT_TYPE input_read(char *input_N){
    char *error_string;
    char *end_pointer;
    INT_TYPE N = 0;
    N = strtol(input_N, &end_pointer, 10);
    if (end_pointer == input_N) {
        error("No digits found, when decoding the input number N of nodes.\n");
    } else if (*end_pointer != '\0') {
        sprintf(error_string, "Invalid digit %c found, when decoding the input number N of nodes.\n", *end_pointer);
        error(error_string);
    }
    if(N>2){
        return N;
    } else {
        error("The number of nodes N must be larger than 2.");
    }
}


//Main
void main(int argc, char **argv){
    
    //Initializing variables
    int i=0, j=0, end_file;
    FILE **fp_array, *fp_single;
    char f_string[20], temp;
    INT_TYPE N = input_read(argv[1]);
    char *data_dir = argv[2];
    INT_TYPE M = 0;
    double c = 0.;
    //Defining a size array with [size_max, size_square_mean]
    INT_TYPE size_array[2];
    size_array[0]=1;
    size_array[1]=N;

    //Initializing options
    #ifdef SEED
    srand(SEED);
    #endif
    #ifndef SEED
    srand(time(0));
    #endif
    char rm_string[4];
    #ifdef WIN_OS
        {sprintf(rm_string, "del");}
    #else
        {sprintf(rm_string, "rm");}
    #endif

    //Starting the code
    printf("\n\nRandom graph analysis with N = %d nodes.\n\n",N);

    //Allocating memory
    nodes=(struct graph *)calloc(N, sizeof(struct graph));
    fp_array=(FILE **)calloc(N_TRAJECTORIES, sizeof(FILE *));

    //Calculating the values for each trajectory
    for(i=0;i<N_TRAJECTORIES;i++){
        //Initializing parameters
        M = 0;
        size_array[0] = 1;
        size_array[1] = N; 
        //Initializing the graph
        initialization(N);
        //Choosing the file name to store the trajectory data 
        sprintf(f_string,"%sgraph%d.dat",data_dir,i);
        //Opening the file to store the trajectory data 
        if((fp_single=fopen(f_string,"w+"))==NULL) error("ERROR: I cannot open (w+) the data files to store the trajectories.");
        //Saving the labels of the data matrix
        fprintf(fp_single,"#Smax       S^2mean       c\n");
        //Saving the trajectory data in the data matrix
        for(c=0.;c<1.;){
            step(N,M,size_array);
            M++;
            c=((double)M/N);
            fprintf(fp_single,"%g   %g  %g\n", (double)size_array[0]/N, (double)(size_array[1]-(INT_TYPE)size_array[0]*size_array[0])/N, c);
            }
        //Closing the file
        fclose(fp_single);
    }

    //Opening the trajectory files in read only mode
    for(i=0;i<N_TRAJECTORIES;i++){
        //Defining the file name
        sprintf(f_string,"%sgraph%d.dat",data_dir,i);
        //Opening and checking open result
        if((fp_array[i]=fopen(f_string,"r"))==NULL){
            printf("\n\nFILE ERROR: %s",f_string); 
            error("ERROR: I cannot open (r) the trajectory files to read them.");
            }
    }

    //Defining the file name of the file where the average quantities willl be saved
    sprintf(f_string,"%sg_mean_%d.dat",data_dir,N);
    //Opening the file with the average quantities in write mode
    if((fp_single=fopen(f_string, "w+"))==NULL) error("ERROR: I cannot open (w+) the file to store the average results.");
    //Saving the labels of the data matrix
    fprintf(fp_single,"Smax_mean,<S^2>_mean,c,dev_Smax_mean,dev_<S^2>_mean\n");
    //Placing the reading buffer to the end of the first line
    for(i=0;i<N_TRAJECTORIES;i++){
        do{
            fscanf(fp_array[i],"%c",&temp);
            }while(temp!='\n');
        }

    //Initializing the average and standard deviations among trajectories
    read_smax.mean=0.; 
    read_ssqmean.mean=0.; 
    read_smax.std_dev=0.; 
    read_ssqmean.std_dev=0.;

    //Computing the average and standard deviations among trajectories
    for(end_file=0; end_file!=EOF; read_smax.mean=0., read_ssqmean.mean=0., read_smax.std_dev=0., read_ssqmean.std_dev=0.){
        //Scrolling the data files
        end_file=fscanf(fp_array[0], "%Lf   %Lf    %Lf", &read_smax.value[0], &read_ssqmean.value[0], &c);
        read_smax.mean+=read_smax.value[0];
        read_ssqmean.mean+=read_ssqmean.value[0];
        for(i=1;i<N_TRAJECTORIES;i++){
            fscanf(fp_array[i], "%Lf   %Lf    %*f", &read_smax.value[i], &read_ssqmean.value[i]);
            read_smax.mean+=read_smax.value[i];
            read_ssqmean.mean+=read_ssqmean.value[i];
        }
        //Computing the average
        read_smax.mean/=(double)N_TRAJECTORIES;
        read_ssqmean.mean/=(double)N_TRAJECTORIES;
        //Evaluating the std of the average
        read_smax.std_dev=std_dev_calc(read_smax);
        read_ssqmean.std_dev=std_dev_calc(read_ssqmean);
        if(end_file!=EOF)
        fprintf(fp_single,"%Lg,%Lg,%Lg,%Lg,%Lg\n", read_smax.mean, read_ssqmean.mean, c, read_smax.std_dev, read_ssqmean.std_dev);
    }

    //Deleting the files with the trajectories (only if instructed to do so)
    for(i=0;i<N_TRAJECTORIES;i++){
        fclose(fp_array[i]);
        #ifdef DELETE
        sprintf(f_string,"%s Data\\graph%d.dat",rm_string, i);
        if(system(f_string)==-1) error("ERROR: I cannot delete the trajectory files");
        #endif
    }
    fclose(fp_single);

}
