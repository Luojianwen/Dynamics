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
    int ret;
    key_t key;
    char *shmadd;
    double *num;
 
    //创建key值
    key = ftok("../", 2015);
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
    num = shmat(shmid, NULL, 0);
    if(shmadd < 0)
    {
        perror("shmat");
        exit(-1);
    }
 
    //读共享内存区数据
    //printf("data = [%s]\n", shmadd);
    while(1){
        printf("data = [%f, %f, %f, %f]\n", num[0], num[1], num[2], num[3]);
    }
 
    /*
    //分离共享内存和当前进程
    //ret = shmdt(shmadd);
    ret = shmdt(num);
    if(ret < 0)
    {
        perror("shmdt");
        exit(1);
    }
    else
    {
        printf("deleted shared-memory\n");
    }
 
    //删除共享内存
    shmctl(shmid, IPC_RMID, NULL);
 
    //system("ipcs -m"); //查看共享内存
    */
    return 0;
}
