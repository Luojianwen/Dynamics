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
    key = ftok("../", 2015);
    if(key == -1)
    {
        perror("ftok");
    }
 
    //创建共享内存
    shmid = shmget(key, BUFSZ, IPC_CREAT|0666);    
    if(shmid < 0)
    {
        perror("shmget");
        exit(-1);
    }
 
    //映射
    //shmadd = shmat(shmid, NULL, 0);
    num = (double*) shmat(shmid, NULL, 0);
    //if(shmadd < 0)
    if(*num < 0)
    {
        perror("shmat");
        _exit(-1);
    }
 
    //拷贝数据至共享内存区
    printf("copy data to shared-memory\n");
    //bzero(shmadd, BUFSZ); // 共享内存清空
    //strcpy(shmadd, "how are you, lh");

    double cnt = 0.1;
    num[0] = 10.0;
    num[1] = 10.1;
    num[2] = 10.2;
    num[3] = 10.3;
    while(1)
    {
        num[0] = num[0] + cnt;
        num[1] = num[1] + cnt;
        num[2] = num[2] + cnt;
        num[3] = num[3] + cnt;
        sleep(1);
    }

    return 0;
}

