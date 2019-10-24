#define MAX_CHILD_PROCESS 8
#define PROJECT_ID 88

int proj_id = PROJECT_ID;
int child_process_fork;
int shared_memory_size[MAX_CHILD_PROCESS];

int ibeg[MAX_CHILD_PROCESS];
int iend[MAX_CHILD_PROCESS];
int jbeg[MAX_CHILD_PROCESS];
int jend[MAX_CHILD_PROCESS];

class SharedMemory
{
private:
  int shmid[MAX_CHILD_PROCESS];
  double *data[MAX_CHILD_PROCESS];
  int nb_child_process_fork;
  // Prevent copying.
  SharedMemory(const SharedMemory&);
  SharedMemory& operator=(const SharedMemory&);

public:
  SharedMemory() {};
  void Init(const int& child_process_fork,
            const int* const shared_memory_size);
  ~SharedMemory();
  int GetShmid(const int& current_fork) const
  {
    return shmid[current_fork];
  }
  double* GetData(const int& current_fork) const
  {
    return data[current_fork];
  }
};

void SharedMemory::Init(const int& nb_fork, const int* const shared_mem_size)
{
  nb_child_process_fork = nb_fork;
  // Set shared memory.
  for (int i = 0;  i < nb_child_process_fork - 1; i++)
    {
      char tmpfile[] = "/tmp/fileXXXXXX";
      int fd;

      fd = mkstemp(tmpfile);

      if (fd == -1)
        cerr << "Module: failed to create unique temporary file" << endl;

      key_t cle = ftok(tmpfile, proj_id);

      shmid[i] = shmget(cle, shared_mem_size[i] * sizeof(double),
                        IPC_CREAT | IPC_EXCL | 0666);

      if (shmid[i] == -1)
        cerr << "Module: failed to create unique shared memory segment"
             << endl;

      data[i] = (double*) shmat(shmid[i], NULL, SHM_R);
    }
}

SharedMemory::~SharedMemory()
{
  for (int i = 0; i < nb_child_process_fork - 1; i++)
    {
      shmdt(data[i]);
      if (shmctl(shmid[i], IPC_RMID, NULL) < 0)
        cout << "can not remove shared memory" << shmid[i] << endl;
    }
}
