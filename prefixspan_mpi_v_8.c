/**
*  v8 改进
*改进点：
*   两个串caaaaab caaaaad前6个字符都是完全相同的，在没有加入剪切的机制的时候对caaaaab caaaaad两个串比较完了之后还
*会对aaaaab aaaaad等串进行比较，本次修改就是加入一个检测的机制，检测多个串的前一个分块的值是完全否相同，如果相同则不再
*对他们进行比较，因为这一步分包含在他们的父序列中运行。
*好处：
*   减少了比较次数，并且减少了重复输出的序列数
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include <math.h>
#include <getopt.h>

#define  SAVA_RESULT_NUM  5000         //保存生成结果时一次分配的空间大小
#define  SAVA_ADD_NUM 1000				     //每次申请的nodeaddress数量
#define  MAX_BUF_SIZE 409600           //设置缓存的长度
#define  MAX_FDATA_NODE_LEN 409600
#define  NODE_NUM 4                    //基因种类一般只有4种
#define  NODE_LEN 10                  //将数据中长度为NODE_NUM的子串作为一个元素
#define  ONECE_S_NUM 200000            //一次生成投影数据库的节点数
#define  ONECE_S_NUM_IN_MIAN  2000		 //保存子串的一次性分配的空间个数
/**
* 转换表  tranTable:
*转换表，将基因用数字表示比如基因a可以表示为tranTable['a'-'a']=tranTable[0]=0
*基因'c'可以表示为tranTable['c'-'a']=tranTable[2]=1
*转换表  genne_table:
*将数字表示的基因转换为字符形式
*/
int tranTable[20]={ 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3};
char genne_table[4]={'a','c','g','t'};

/**
* 用于保存文件中的数据
*一条基因存储于一个文件中，将一条基因分散保存于多个基因组中
*处理时还是当作一条基因处理
*/
typedef struct fdata
{
	char c[MAX_FDATA_NODE_LEN];   //长度为MAX_FDATA_NODE_LEN的基因组
	int len;                      //此组基因中基因的条数
	int line;                     //此基因组的编号
	struct fdata *next;           //下一个基因组的地址
} fdata;

/// 用来保存子串地址
typedef struct  nodeaddress
{
	int distancetostrart;         //相对于此基因组开始处一偏移量
	int sd;                       //此条子串的起始处偏移量
	int ed;                       //此条子串的结尾处偏移量
	int parent;                  //父串的值
	struct fdata *line;           //节点所处的基因组的地址
	struct nodeaddress *next;     //此条子串下一处出现的位置
}nodeaddress;
/**
* 此结构体用来保存一组用于保存子串地址
*一次性分配SAVA_ADD_NUM个地址结构体
*/
typedef struct sava_nodeaddress
{
	struct nodeaddress add[SAVA_ADD_NUM];  //保存SAVA_ADD_NUM个地址
	struct sava_nodeaddress *next;         //下一个
	struct nodeaddress *usenow;            //在将地址写入时，当前已经使用到的地址
}sava_nodeaddress;

/// 用于保存满足支持度的子串
typedef struct cutednode
{
	int nodevalue;                         //子串的值
	int num;                               //子串的支持数
	int len;                               //子串的长度
	struct nodeaddress *firstaddresss;     //子串第一次出现的位置，此是一个链表，保存全部出现的位置
	struct nodeaddress *now;               //在向地址链表写入地址时，当前用到的位置
	struct cutednode *next;                //下一个子串
}cutednode;

/// 在初始化时用于保存头次扫描所有满足条件的子串
typedef struct sava_cutednode
{
	struct cutednode cutname[ONECE_S_NUM];
	struct sava_cutednode *next;               //下一个ONECE_S_NUM条子串的位置
	struct cutednode *nowuse;                  //当进行运行时下在运行的子串
}sava_cutednode;

typedef struct sava_cutednode_in_main
{
	struct cutednode cutname[ONECE_S_NUM_IN_MIAN];   //保存非初始子串的ONECE_S_NUM_IN_MIAN条子串
	struct sava_cutednode_in_main *next;             //下一个ONECE_S_NUM_IN_MIAN条子串
	struct cutednode *nowuse;                        //现在正在使用的非初始子串
}sava_cutednode_in_main;

///用来保存节点
typedef struct node
{
	int nodevalue;                                //四叉树的节点值
	int num;                                      //以此节点结束的串的支持数
	int trun;                                     //运行的轮数
	int sd;                                       //保存此串的头元素出现的位置
	int parent;                                  //此数据用来保存此子串的父串的值，如果这此串所有出现位置的父节点相同则保存父节点的值，否则值为-1
	int ed;										 //保存最后一次将子串的
	struct nodeaddress *fistadd;                  //保存此节点所有出现的位置，此指针指向第一个位置
	struct nodeaddress *addnow;					          //保存最近一将次出现的位置
    struct cutednode *sava;                       //此节点对应的子串的地址
	struct node *next;
}node;

///用于保存当前运行的子串和子串的子串运行状态
typedef struct mainrun_sava_onece_node
{
	struct cutednode *first;                            //保存此阶段的非初始子串的此组中的起始位置
 	struct cutednode *nowuse;                           //现在正在运行的子串的地址
	struct mainrun_sava_onece_node *next;               //下一阶段的子串
	struct   mainrun_sava_onece_node  *up;              //上一阶段的子串的位置
	int len;
}mainrun_sava_onece_node;

void readfile (char *path);    //读取文件中的数据保存于结构体fdata中
void sava_node(node *temp,int len,sava_cutednode_in_main **cutednode_f,int *all_cutednode_num,int *sava_cutednode_num,cutednode **recover_table);
void malloc_sava_add_node();
void init_node(int myid,int numprocs);  //初始化数据取出所有长度为NODE_LEN的所有子串
void create_node(cutednode *nodetemp,int trun,int num);
void sava_cutednode_address(sava_nodeaddress **sava_addreass_first,sava_nodeaddress **sava_addreass_now);
void get_shadow_data(cutednode *cutednode_temp,int now_use_num );
void pop_node_to_node(int a, int d,fdata *fdata_temp,int s,int len,int sd,int parent);
void get_node(nodeaddress *addtemp,int len,int *mainrun_cutednode_allnum,int *mainrun_cutednode_num,cutednode **recover_table);
void mainrun(int numprocs,int myid,char *processor_name,char *savapath);
int main(int argc,char *argv[]);
int chartoint(char str[]);

int min_num,min_node_len_in_this;   //支持度,最小长度
sava_cutednode *sava_cutednode_first,*sava_cutednode_now;
sava_cutednode_in_main *ssava_cutednode_first,*ssava_cutednode_now;

int sava_cutednode_num=0,sava_add_num=0;
int fsava_nodeaddress_num=0,ssava_nodeaddress_num=0;
int fsava_nodeaddress_allnum=0,ssava_nodeaddress_allnum=0;
sava_nodeaddress *fsava_nodeaddress_first=NULL,*fsava_nodeaddress_now,*ssava_nodeaddress_first=NULL,*ssava_nodeaddress_now;
cutednode *cutednode_now,*cutednode_first;

nodeaddress *fnodeaddress_first,*fnodeaddress_now;
int size,trun=1,init_cutednode_usenum=0;
node *node_first=NULL,*last_node_first=NULL,*now;
fdata *fdata_first,*t;
int cut_value;
int main(int argc,char *argv[])
{

	FILE *fp;
	int namelen,myid,numprocs,allnum=0,minlen=0,minnum=0,next_option;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	double startwtime,time1,time2,maxtime,this_node_use_time;
	char *path=NULL,*savapath=NULL,*len=NULL,*num=NULL;
	const char* const short_options="hp:s:l:n:";
	const struct option long_options[]={{"path",1,NULL,'p'},{"savapath",1,NULL,'s'},{"len",1,NULL,'l'},{"num",1,NULL,'n'},{"help",1,NULL,'h'}};
	 cut_value=(int)pow(4.0,NODE_LEN)-1	;
	do{
		next_option=getopt_long(argc,argv,short_options,long_options,NULL);
		switch(next_option)
		{
			case 'h':
				printf("  -p  --path:\n        数据文件路径\n"
				       "  -s  --savapath:\n        "
				       "结果保存文件的路径\n "
				       " -l  --len:\n        最小长度\n "
				       " -n  --num:\n        支持度\n "
				       " -h  --help:\n        帮助\n");
				return 0;
				break;
			case 'p':
				path=optarg;
				break;
			case 's':
				savapath=optarg;
				break;
			case 'l':
				len=optarg;
				break;
			case 'n':
				num=optarg;
				break;
		}
	}while(next_option!=-1);
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);
    startwtime = MPI_Wtime();

    if(myid==0){
		if(len==NULL || num==NULL || path==NULL || savapath==NULL){
			printf("error options\n");
			MPI_Abort(MPI_COMM_WORLD,1);
			return 0;
		}
	}

	min_node_len_in_this=chartoint(len);
	min_num=chartoint(num);

	if(myid==0){
		if(min_num<1 || min_node_len_in_this <10){//此处暂时设置最小长度为10
			printf("error options\n");
			MPI_Abort(MPI_COMM_WORLD,1);
			return 0;
		}
	}

	readfile(path);
	init_node(myid,numprocs);
	time1=MPI_Wtime();
	printf("%d %s \tinit_node  and readfile use time:%f \n",myid,processor_name,time1-startwtime);

	mainrun(numprocs,myid,processor_name,savapath);
	time2=MPI_Wtime();
	printf("%d %s \tmainrun use time:%f  \tall time:%f\n",myid,processor_name,time2-time1,time2-startwtime);

	MPI_Reduce(&init_cutednode_usenum,&allnum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	if(myid==0){
		printf("all line num %d \tmax use time: %f \n",allnum,MPI_Wtime()-startwtime);
    }
	MPI_Finalize();
	return 0;
}

///将字符串转换为int类型
int chartoint(char str[])
{
	int len=0,i,num=0;
	len=strlen(str);
	for(i=0;i<len;i++){
		if(str[i]<48 || str[i]>57){
			printf("error  option\n");
			MPI_Abort(MPI_COMM_WORLD,2);
			return 0;
		}
		num=10*num+(str[i]-'0');
	}
	return num;
}

/**
* 将支持数大于支持度的元素保存
*temp:四叉树的尾部
*recover_table:回收表，在非初始化子串运行时如果子串不满足条件则为它分配
*的空间放入回收表中以供其它子串使用
*/
void sava_node(node *temp,int len,sava_cutednode_in_main **cutednode_f,int *all_cutednode_num,int *sava_cutednode_num,cutednode **recover_table)
{
	cutednode *cutednode_temp;
	node *node_temp;
	sava_cutednode_in_main *sava_cutednode_temp;
	int i;
	i=len-1;
	node_temp=temp;
	if((*recover_table)!=NULL){//判断回收表中是否有空间供使用
		cutednode_temp=(*recover_table);
		(*recover_table)=(*recover_table)->next;
	}else{
		if((*sava_cutednode_num)%ONECE_S_NUM_IN_MIAN==0){
			if((*sava_cutednode_num)>=(*all_cutednode_num)){
				if((sava_cutednode_temp=(sava_cutednode_in_main *)malloc(sizeof(sava_cutednode_in_main)))==NULL){
					printf("error in malloc sava_cutednode \n");
					MPI_Abort(MPI_COMM_WORLD,21);
					exit(21);

				}
				sava_cutednode_temp->next=NULL;
				sava_cutednode_temp->nowuse=sava_cutednode_temp->cutname;
				if((*cutednode_f)==NULL){
					(*cutednode_f)=sava_cutednode_temp;
				}else{
					ssava_cutednode_now->next=sava_cutednode_temp;
				}
				ssava_cutednode_now=sava_cutednode_temp;
				(*all_cutednode_num)+=ONECE_S_NUM_IN_MIAN;
			}else{
				if((*sava_cutednode_num)==0) ssava_cutednode_now=(*cutednode_f);
				else ssava_cutednode_now=ssava_cutednode_now->next;
				ssava_cutednode_now->nowuse=ssava_cutednode_now->cutname;
			}


		}
		cutednode_temp=ssava_cutednode_now->nowuse;
		ssava_cutednode_now->nowuse++;
	}

	(*sava_cutednode_num)++;
	cutednode_temp->len=len;
	cutednode_temp->firstaddresss=node_temp->fistadd;
	cutednode_temp->num=node_temp->num;
	cutednode_temp->next=NULL;
	if(cutednode_first==NULL){
		cutednode_first=cutednode_temp;
	}else{
		cutednode_now->next=cutednode_temp;

	}
	cutednode_now=cutednode_temp;
	cutednode_temp->nodevalue=node_temp->nodevalue;
}

/// 程序初始化时取出所有可行节点                                          */
void init_node(int myid,int numprocs)
{
	fdata *fdata_temp;
	int i=0,fdata_len,d=0,t,sava_cutednode_num_trun=0,this_sequnce_num=0,last_sequnce_num=-1;
	char *thisnode;
	node *add,*temp;
	cutednode *cutednode_temp=NULL;
	sava_cutednode *sava_cutednode_temp;
	fdata_temp=fdata_first;

	/****************创建hash表****************/
	if((node_first=(node *)malloc(sizeof(node)*(cut_value+1)))==NULL){
		printf("error in malloc node in init_node\n");
		MPI_Abort(MPI_COMM_WORLD,1);
		exit(0);
	}
	/****************创建完成*********************/

	while(fdata_temp!=NULL)
	{
		i=0;
		fdata_len=fdata_temp->len;
		thisnode=fdata_temp->c;
		while(i<fdata_len)
		{
			/**当此数存在**/
			if(*thisnode>96  && *thisnode<117) {
					d++;
					t=tranTable[*thisnode-'a'];
					this_sequnce_num<<=2;//向左移动两位因为基因只有4种
					this_sequnce_num+=t;
					this_sequnce_num=this_sequnce_num&cut_value;

					/**pop_node_to_node_first   start
					* 初始化时，将长度为r为的子串保存，如果已经存在则将支持数加1
					* a:要保存的子串
					* j:在数组a中，子串的起始地址
					* d:子串尾部离文件起始位置的距离
					*/
					if(d>=NODE_LEN){
								temp=&node_first[this_sequnce_num];
								if(temp->trun==trun){
                  if(temp->parent!=last_sequnce_num) temp->parent=-1;
									if(d-(temp->ed)>=NODE_LEN){
												temp->ed=d;
												temp->sd=d-NODE_LEN;
												temp->num++;
									}
								}else{
									temp->nodevalue=this_sequnce_num;
									temp->ed=d;
									temp->sd=d-NODE_LEN;
									temp->next=NULL;
									temp->num=1;
									temp->trun=trun;
                                    temp->fistadd=NULL;
                                    temp->parent=last_sequnce_num;
                                    if(last_node_first!=NULL){
                                        now->next=temp;
                                    }else{
                                        last_node_first=temp;
                                    }
                                    now=temp;
							 }
							 last_sequnce_num=this_sequnce_num;
					}
					/**pop_node_to_node_first  end*/
			}
			i++;
			if(i!=fdata_len) thisnode++;
		}
		fdata_temp=fdata_temp->next;
	}
	sava_cutednode_first=NULL;
	sava_cutednode_now=NULL;

	temp=last_node_first;
	while(temp!=NULL)
	{
		//当支持数大于支持度
		if(temp->num>=min_num ){
			/**sava_node  函数 start*/
			if(sava_cutednode_num_trun==myid){
				i=NODE_LEN-1;
				add=temp;
				if(init_cutednode_usenum%ONECE_S_NUM==0){
					if((sava_cutednode_temp=(sava_cutednode *)malloc(sizeof(sava_cutednode)))==NULL){
							printf("error in malloc sava_cutednode \n");
							MPI_Abort(MPI_COMM_WORLD,21);
							exit(21);
						}
						sava_cutednode_temp->next=NULL;
						sava_cutednode_temp->nowuse=sava_cutednode_temp->cutname;
						if(sava_cutednode_first==NULL){
							sava_cutednode_first=sava_cutednode_temp;
						}else{
							sava_cutednode_now->next=sava_cutednode_temp;
						}
						sava_cutednode_now=sava_cutednode_temp;
				}
				cutednode_temp=sava_cutednode_now->nowuse;
				sava_cutednode_now->nowuse++;

				init_cutednode_usenum++;
				cutednode_temp->len=NODE_LEN;
				cutednode_temp->firstaddresss=add->fistadd;
				cutednode_temp->num=add->num;
				cutednode_temp->next=NULL;
				if(cutednode_first==NULL){
					cutednode_first=cutednode_temp;
				}else{
					cutednode_now->next=cutednode_temp;
				}
				cutednode_now=cutednode_temp;

				cutednode_temp->nodevalue=temp->nodevalue;
				/**over*/
			}
			sava_cutednode_num_trun++;
			if(sava_cutednode_num_trun==numprocs) sava_cutednode_num_trun=0;
		}
		temp=temp->next;
	}
}

///从文件中获取数据
void readfile(char *path)
{
	fdata *fdata_temp,*fdata_now;
	char *buf,*s;
	int len,i;
	FILE *datafile;

	if((datafile=fopen(path,"r"))==NULL)
	{
		printf("error in open file %s\n",path);
		MPI_Abort(MPI_COMM_WORLD,13);
		exit(13);
	}

	if((buf=(char *)malloc(MAX_BUF_SIZE*sizeof(char)))==NULL){
		printf("error in malloc buf \n");
		MPI_Abort(MPI_COMM_WORLD,4);
		exit(4);
	}
	if((s=(char *)malloc(MAX_BUF_SIZE*sizeof(char)))==NULL){
		printf("error in malloc s \n");
		MPI_Abort(MPI_COMM_WORLD,5);
		exit(5);
	}
	setvbuf(datafile,buf,_IOFBF,MAX_BUF_SIZE);
	memset(s,'\0',MAX_BUF_SIZE);
	i=0;
	while (fgets(s,MAX_BUF_SIZE,datafile))
	{
		len=strlen(s);
		if(len<15) continue;//因为长度小于15有可能是文件中的注释

		if((fdata_temp=(fdata *)malloc(sizeof(fdata)))==NULL)
		{
			printf("error in read file\n");
			MPI_Abort(MPI_COMM_WORLD,6);
			exit(6);
		}
		strncpy(fdata_temp->c,s,MAX_BUF_SIZE);
		fdata_temp->len=len;
		fdata_temp->line=i;
		fdata_temp->next=NULL;

		if(fdata_first==NULL){
			fdata_first=fdata_temp;
		}
		else{
			fdata_now->next=fdata_temp;
		}
		fdata_now=fdata_temp;
		memset(s,'\0',MAX_BUF_SIZE);
		i++;
	}
	free(buf);
	free(s);
	fclose(datafile);
}

/**
*如果在第一次扫描整个数据生成所有满足第一次要求的子串时就生成此子串的投影
* 数据库这将花费大量时间，所以在第一次扫描 的时候并不生成投影数据 库，而
* 一在第一次运行此串时生成但是如果每次扫描文件只为一个串生成投影数据库，
*这将花费大量时间，所以每次为ONCE_S_NUM个串生成子投影数据库
*/
void create_node(cutednode *nodetemp,int trun,int num)
{
	node *temp;
	cutednode *cutednode_this_temp;
	int i=0,j;
	cutednode_this_temp=nodetemp;
	while(num>0 && cutednode_this_temp!=NULL)
	{
		temp=&node_first[cutednode_this_temp->nodevalue];
		temp->trun=trun;
		temp->num=0;

		temp->sava=cutednode_this_temp;

		cutednode_this_temp=cutednode_this_temp->next;
		num--;
  }
}
void sava_cutednode_address(sava_nodeaddress **sava_addreass_first,sava_nodeaddress **sava_addreass_now)
{
	sava_nodeaddress *sava_nodeaddress_temp;
	if((sava_nodeaddress_temp=(sava_nodeaddress *)malloc(sizeof(sava_nodeaddress)))==NULL){
		printf("error in sava_cutednode_address malloc sava_nodeaddress_temp\n");
		MPI_Abort(MPI_COMM_WORLD,7);
		exit(7);
	}
	sava_nodeaddress_temp->next=NULL;

	if((*sava_addreass_first)==NULL){
		(*sava_addreass_first)=sava_nodeaddress_temp;
	}else{
		(*sava_addreass_now)->next=sava_nodeaddress_temp;
	}
	(*sava_addreass_now)=sava_nodeaddress_temp;

}

/**
* 获取投影数据库
* nodetemp:子串的值
* addtemp: 用来保存此子投影数据库
*/
void get_shadow_data(cutednode *cutednode_temp,int now_use_num )
{
	node *temp,*node_temp[NODE_LEN];
	fdata *fdata_temp;
	int i,j=0,mark=0,fdata_len,d=0,n=0,t,this_sequnce_num=0,last_sequnce_num=-1;

	char *thisnode,nodevalue[NODE_LEN];
	cutednode *nodetemp;
	nodeaddress *addtemp1;

	nodetemp=cutednode_temp;
	trun++;
	create_node(nodetemp,trun,now_use_num);

	fsava_nodeaddress_num=0;

	fdata_temp=fdata_first;
	while(fdata_temp!=NULL)
	{
		i=0;
		fdata_len=fdata_temp->len;
		thisnode=fdata_temp->c;
		while(i<fdata_len)
		{
			if(*thisnode!='\n' && *thisnode!=' ' && *thisnode!='\0') {
					d++;
					t=tranTable[*thisnode-'a'];
					this_sequnce_num<<=2;//向左移动两位因为基因只有4种
					this_sequnce_num+=t;
					this_sequnce_num=this_sequnce_num&cut_value;

					if(d>=NODE_LEN){
						/****************************************/
                            temp=&node_first[this_sequnce_num];
                            if(temp->trun==trun){
                                if(fsava_nodeaddress_num%SAVA_ADD_NUM==0){
                                    if(fsava_nodeaddress_num>=fsava_nodeaddress_allnum){
                                        sava_cutednode_address(&fsava_nodeaddress_first,&fsava_nodeaddress_now);
                                        fsava_nodeaddress_allnum+=SAVA_ADD_NUM;
                                    }else{
                                        if(fsava_nodeaddress_num==0) fsava_nodeaddress_now=fsava_nodeaddress_first;
                                        else fsava_nodeaddress_now=fsava_nodeaddress_now->next;
                                    }
                                    fsava_nodeaddress_now->usenow=fsava_nodeaddress_now->add;
                                }
                                addtemp1=fsava_nodeaddress_now->usenow;
                                fsava_nodeaddress_now->usenow++;
                                fsava_nodeaddress_num++;

                                addtemp1->line=fdata_temp;
                                addtemp1->distancetostrart=i;
                                addtemp1->ed=d;
                                addtemp1->sd=d-NODE_LEN;
                                addtemp1->next=NULL;
                                addtemp1->parent=last_sequnce_num;

                                if(temp->sava->firstaddresss==NULL) {
                                    temp->sava->firstaddresss=addtemp1;
                                }else{
                                    temp->sava->now->next=addtemp1;
                                }
                                temp->sava->now=addtemp1;
                            }
                            last_sequnce_num=this_sequnce_num;
							/******************************************************************/
					}
			}
			i++;
			if(i!=fdata_len) thisnode++;
		}
		fdata_temp=fdata_temp->next;
	}
}

/**
* 将子串加入四叉树中
* a:要为入的子串的值
* d:离整个文件开始处的距离
* fdata_temp :此子串所在的行
* s:子串最后一个值离行开始处的距离
*/
void pop_node_to_node(int a, int d,fdata *fdata_temp,int s,int len,int sd,int parent)
{

	node *temp,*add;
	int i=0,mark=0,t;
	nodeaddress *addtemp;
	sava_nodeaddress *tt;
	tt=ssava_nodeaddress_first;

	temp=&node_first[a];

	if(ssava_nodeaddress_num%SAVA_ADD_NUM==0){
		 if(ssava_nodeaddress_num>=ssava_nodeaddress_allnum ){
				sava_cutednode_address(&ssava_nodeaddress_first,&ssava_nodeaddress_now);
				ssava_nodeaddress_allnum+=SAVA_ADD_NUM;
		 }else{
				if(ssava_nodeaddress_num!=0)  ssava_nodeaddress_now=ssava_nodeaddress_now->next;
                else ssava_nodeaddress_now=ssava_nodeaddress_first;
        }
        ssava_nodeaddress_now->usenow=ssava_nodeaddress_now->add;
	}

	addtemp=ssava_nodeaddress_now->usenow;
	ssava_nodeaddress_now->usenow++;
	ssava_nodeaddress_num++;

	addtemp->distancetostrart=s;
	addtemp->line=fdata_temp;
	addtemp->ed=d;
	addtemp->sd=sd;
	addtemp->next=NULL;
	addtemp->parent=parent;

	if(temp->trun==trun){
		 temp->addnow->next=addtemp;
		 temp->addnow=addtemp;
         if(temp->parent!=parent)  temp->parent=-1;
		 if(sd>=(temp->ed)){
					temp->ed=d;
					temp->sd=sd;
					temp->num++;
        }
	 }else{
				temp->ed=d;
				temp->sd=sd;
				temp->num=1;
				temp->parent=parent;
				temp->trun=trun;
				temp->next=NULL;
				temp->fistadd=addtemp;
				temp->addnow=addtemp;
				if(last_node_first!=NULL){
										now->next=temp;
				}else{
										last_node_first=temp;
				}
				now=temp;
	 }
}
/**
* 根据投影数据库生成下一次运行的所有元素
* addtemp:上次生成的投影数据库
*/
void get_node(nodeaddress *addtemp,int len,int *mainrun_cutednode_allnum,int *mainrun_cutednode_num,cutednode **recover_table)
{
	nodeaddress *add;
	fdata *fdata_temp;
	char *thisnode;
	int i,fdata_len,d,j,sd,this_sequnce_num=0,t;
	last_node_first=NULL;
	node *temp;
	add=addtemp;
	trun++;

	while(add!=NULL)
	{
		j=0;
		this_sequnce_num=0;
		i=add->distancetostrart;
		fdata_temp=add->line;
		fdata_len=fdata_temp->len;
		d=add->ed;
		sd=add->sd;
		thisnode=fdata_temp->c+i;
		while(j<len)
		{
			i++;
			if(i!=fdata_len){
				thisnode++;
			}else{
				fdata_temp=fdata_temp->next;
				if(fdata_temp==NULL) break;
				fdata_len=fdata_temp->len;
				thisnode=fdata_temp->c;
				i=0;
			}
			t=*thisnode-'a';
			if(t>-1 && t<20) {
				j++;
				d++;
				this_sequnce_num=this_sequnce_num*4+tranTable[t];//向左移动两位因为基因只有4种
			}
		}
		if(fdata_temp==NULL) break;
		pop_node_to_node(this_sequnce_num,d,fdata_temp,i,len,sd,add->parent);
		add=add->next;
	}
	cutednode_first=NULL;
	temp=last_node_first;

	while(temp!=NULL){
		if(temp->num>=min_num && temp->parent==-1){
				sava_node(temp,len,&ssava_cutednode_first,mainrun_cutednode_allnum,mainrun_cutednode_num,recover_table);//
		}
		temp->addnow=NULL;
		temp->fistadd=NULL;
		temp=temp->next;
	}

}

/// 主运行函数                                                            */
void mainrun(int numprocs,int myid,char *processor_name,char *savapath)
{
	FILE *sava;
	cutednode *nodetemp,*nodetemp1,*nodetemp2,*recover_table=NULL,*recover_table_now;
	int alllen=0,run_line=0,i,mainrun_cutednode_allnum=0,mainrun_cutednode_num=0,mark=0,cutednode_value,result_all_num=SAVA_RESULT_NUM;
	mainrun_sava_onece_node *mainrun_sava_onece_node_first=NULL,*mainrun_sava_onece_node_now=NULL,*mainrun_sava_onece_node_temp;
	char *theresult,*result_temp;
	sava=fopen(savapath,"w");
	nodetemp=cutednode_first;
	double s_start,r_start,s_all=0,r_all=0,g_all=0,g_start;

	if((theresult=(char *)malloc(SAVA_RESULT_NUM*sizeof(char)))==NULL){
		printf("error in malloc theresult_s\n");
		MPI_Abort(MPI_COMM_WORLD,25);
		exit(25);
	}
	while(nodetemp!=NULL)
	{
		s_start=MPI_Wtime();
		get_shadow_data(nodetemp,ONECE_S_NUM);//由于在第一次分离出所有支持数满足支持度的子串时没有生成投影数据库，所以当作特殊处理
		s_all+=MPI_Wtime()-s_start;
		printf("%d %s \thas run line %d  \tall line num is %d*********\n",myid,processor_name,run_line,init_cutednode_usenum);
		ssava_nodeaddress_num=0;
		nodetemp1=nodetemp;
		// 这里这样处理必须满足的一个条件是，最小子串长度（MIN_NODE_LEN_IN_THIS）大于NODE_LEN
        r_start=MPI_Wtime();
		while(1)
		{
			while(nodetemp1!=NULL)
			{
				if(mainrun_sava_onece_node_now==NULL)
				{
					alllen=0;
					run_line++;//运行的行数
				}
				if((1+alllen)>=result_all_num){
					result_all_num+=SAVA_RESULT_NUM;
					if((result_temp=realloc(theresult,result_all_num))==NULL){
						printf("error in realloc result_temp\n");
						MPI_Abort(MPI_COMM_WORLD,26);
						exit(26);
					}
					theresult=result_temp;
				}

				i=alllen+nodetemp1->len-1;
				cutednode_value=nodetemp1->nodevalue;
				while(i>=alllen){
					theresult[i]=genne_table[cutednode_value&3];
					cutednode_value>>=2;
					i--;
				}
				alllen+=nodetemp1->len;
				theresult[alllen]='\0';

                //生成子投影数据库
                g_start=MPI_Wtime();
				sava_cutednode_num=0;
				if(alllen>=min_node_len_in_this){
					 get_node(nodetemp1->firstaddresss,1,&mainrun_cutednode_allnum,&mainrun_cutednode_num,&recover_table);
				}else if(alllen+NODE_LEN<=min_node_len_in_this){
					 get_node(nodetemp1->firstaddresss,NODE_LEN,&mainrun_cutednode_allnum,&mainrun_cutednode_num,&recover_table);
				}else{
						get_node(nodetemp1->firstaddresss,min_node_len_in_this-alllen,&mainrun_cutednode_allnum,&mainrun_cutednode_num,&recover_table);
				}
				g_all+=(MPI_Wtime()-g_start);


				//如果没有再满足的子串如下处理
				if(cutednode_first==NULL){

					if(alllen>=min_node_len_in_this)  fprintf(sava,"%s : %d\n",theresult,nodetemp1->num);
					//如果此串没有满足支持度的子串就将结写入文件，并将些节点释放
					alllen-=nodetemp1->len;
					theresult[alllen]='\0';

					nodetemp1->firstaddresss=NULL;

					// 将空间回收到回收表                                                   */
					if(mainrun_sava_onece_node_now!=NULL){
						mainrun_cutednode_num--;//已经使用的cutednode空间数量减1
						if(recover_table==NULL){
							recover_table=nodetemp1;
						}else{
							recover_table_now->next=nodetemp1;
						}
						recover_table_now=nodetemp1;//这个空间内容已经无效所以将其回收以供其它子串使用
					}

					/**
					* 当一个sava_cutednode中的cutednode运行完成可能出现三种情况
					* 1.此层次所有元素已经运行完成，则返回上一层，运行上一层没有运行完成的子串
					* 2.此sava_cutednode还有子sava_cutednode，则表示此层还没有运行完成，则进入子sava_cutednode继续运行
					* 3.此sava_cutednode还没有运行完成，则继续在此sava_cutednode中运行
					*/
					if(nodetemp1->next==NULL || mainrun_sava_onece_node_now==NULL){
						break;
					}else{
						nodetemp2=nodetemp1->next;
						nodetemp1->next=NULL;
						nodetemp1=nodetemp2;

						mainrun_sava_onece_node_now->nowuse=nodetemp1;
					}
				}else{

					if(mainrun_sava_onece_node_now!=NULL){
						if(mainrun_sava_onece_node_now->next==NULL){
							if((mainrun_sava_onece_node_temp=(mainrun_sava_onece_node *)malloc(sizeof(mainrun_sava_onece_node)))==NULL){
								printf("error in malloc mainrun_sava_onece_node_temp \n");
								MPI_Abort(MPI_COMM_WORLD,30);
								exit(30);
							}
							mainrun_sava_onece_node_temp->next=NULL;
							mainrun_sava_onece_node_now->next=mainrun_sava_onece_node_temp;
							mainrun_sava_onece_node_temp->up=mainrun_sava_onece_node_now;
							mainrun_sava_onece_node_now=mainrun_sava_onece_node_temp;
						}else{
							mainrun_sava_onece_node_now=mainrun_sava_onece_node_now->next;
						}
					}else{
						if(mainrun_sava_onece_node_first==NULL){
							if((mainrun_sava_onece_node_first=(mainrun_sava_onece_node *)malloc(sizeof(mainrun_sava_onece_node)))==NULL){
								printf("error in malloc mainrun_sava_onece_node_temp \n");
								MPI_Abort(MPI_COMM_WORLD,30);
								exit(30);
							}
							mainrun_sava_onece_node_first->next=NULL;
							mainrun_sava_onece_node_first->up=NULL;

						}
						mainrun_sava_onece_node_now=mainrun_sava_onece_node_first;
					}
					mainrun_sava_onece_node_now->first=cutednode_first;
					mainrun_sava_onece_node_now->nowuse=cutednode_first;
					mainrun_sava_onece_node_now->len=alllen;
					nodetemp1=cutednode_first;
					continue;
				}
			}

			if(mainrun_sava_onece_node_now!=NULL){
				//表示还有非初始子串没有运行完成                                       */
				mainrun_sava_onece_node_now=mainrun_sava_onece_node_now->up;
				if(mainrun_sava_onece_node_now!=NULL){

					alllen=mainrun_sava_onece_node_now->len;
					theresult[alllen]='\0';

					mainrun_sava_onece_node_now->nowuse=mainrun_sava_onece_node_now->nowuse->next;
					nodetemp1=mainrun_sava_onece_node_now->nowuse;
				}else{
					nodetemp1=NULL;
				}
			}else{
				/**
				* 当一个初始子串运行完成
				*可能出现两种情况：
				*1.此ONECE_S_NUM个初始子串还没有运行完成，则继续运行下一个初始子串
				*2.运行完成则运行下一次初始投影数据库生成
				*如果为2则将sava_cutednode_now赋值为ssava_cutednode_first，并退出此次循环
				*/
				ssava_nodeaddress_num=0;
				mainrun_cutednode_num=0;
				recover_table=NULL;
				ssava_cutednode_now=ssava_cutednode_first;
				if(ssava_cutednode_first!=NULL) ssava_cutednode_now->nowuse=ssava_cutednode_first->cutname;
				nodetemp=nodetemp->next;
				if(nodetemp==NULL || nodetemp->firstaddresss==NULL){
						break;
				}
				nodetemp1=nodetemp;
			}
		}
		r_all+=MPI_Wtime()-r_start;
		fflush(sava);
	}
	printf("%d %s \trun: %f    \tget_shadow_data:%f  \tget_node:%f \t\n",myid,processor_name,r_all,s_all,g_all);
	fclose(sava);
}
