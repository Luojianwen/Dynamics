#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#define BUFSZ 32

int main(int argc, char *argv[])
{
    int shmid;
    key_t key;
    //char *shmadd; // for char variable
    double *num;
 
    //创建key值
    key = ftok("~/", 2015); // 1st parameter is file path, should be aligned with writing procedure.
    if(key == -1)
    {
        perror("ftok");
    }
    
    //system("ipcs -m"); //查看共享内存
    
    //打开共享内存
    shmid = shmget(key, BUFSZ, IPC_CREAT|0666);
    if(shmid < 0)
    {
        perror("shmget");
        exit(-1);
    }
 
    //映射
    //shmadd = shmat(shmid, NULL, 0);
    num = (double*) shmat(shmid, NULL, 0);
    if(*num < 0)
    {
        perror("shmat");
        exit(-1);
    }
 
    //读共享内存区数据
    //printf("data = [%s]\n", shmadd);
    while(1){
        printf("data = [%f, %f, %f, %f]\n", num[0], num[1], num[2], num[3]);
        sleep(1);
    }
    
    return 0;
}
