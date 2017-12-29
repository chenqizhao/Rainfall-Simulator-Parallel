#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

#define NUM_THREADS 8  // define # of threads used

int N; //the input grid size
int t = 1;
int rain_time;
double absorb_rate;
bool wet = true;

std::vector<std::vector<int>> land; //landscape
std::vector<std::vector<double> > trap_water;
std::vector<std::vector<std::string> > lowest_neighbours; //lowest neighbors
std::vector<std::vector<double> > absorbed; //absorbed rainfall
std::vector<std::vector<double> > flux; // trickels flow to neighbors

void init() {
  land = std::vector<std::vector<int> > (N, std::vector<int>(N, 0));
  trap_water = std::vector<std::vector<double> > (N, std::vector<double>(N, 0));
  lowest_neighbours = std::vector<std::vector<std::string> > (N, std::vector<std::string>(N, ""));
  absorbed = std::vector<std::vector<double> > (N, std::vector<double>(N, 0));
  flux = std::vector<std::vector<double> > (N, std::vector<double>(N, 0));
}

double calc_time(struct timespec start, struct timespec end) {
  double start_sec = (double)start.tv_sec*1000000000.0 + (double)start.tv_nsec;
  double end_sec = (double)end.tv_sec*1000000000.0 + (double)end.tv_nsec;

  if (end_sec < start_sec) {
    return 0;
  } else {
    return end_sec - start_sec;
  }
}

bool dry() {  // 这个函数可以pthread
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (trap_water[i][j] != 0) {
	return false;
      }
    }
  }
  return true;
}

void* check_dry(void* arg) {  // 这个函数可以pthread
  int id = *(int *)arg;
  int start = id * (N / NUM_THREADS);
  int end = (id + 1) * (N / NUM_THREADS);
  for (int i = start; i < end; i++) {
    for (int j = 0; j < N; j++) {
      if (trap_water[i][j] != 0) {
        wet = true;
        return NULL;
      }
    }
  }
  wet = false;
}

std::vector<int> get_neighbors_heights(int row, int col) {
  int u = 0;
  int d = 0;
  int r = 0;
  int l = 0;

  if (row == 0) {
    u = -1;
    d = land[row + 1][col];
  } else if (row == N - 1) {
    u = land[row - 1][col];
    d = -1;
  } else {
    u = land[row - 1][col];
    d = land[row + 1][col];
  }

  if (col == 0) {
    l = -1;
    r = land[row][col + 1];
  } else if (col == N - 1) {
    l = land[row][col - 1];
    r = -1;
  } else {
    l = land[row][col - 1];
    r = land[row][col + 1];
  }
    
  std::vector<int> ret;
  ret.push_back(u);
  ret.push_back(d);
  ret.push_back(l);
  ret.push_back(r);
  return ret;
}

std::vector<int> find_direction(char ch) {
  int dx = 0;
  int dy = 0;
  if (ch == 'u') {
    dx = -1;
    dy = 0;
  } else if (ch == 'd') {
    dx = 1;
    dy = 0;
  } else if (ch == 'l') {
    dx = 0;
    dy = -1;
  } else if (ch == 'r') {
    dx = 0;
    dy = 1;
  }
  std::vector<int> ret;
  ret.push_back(dx);
  ret.push_back(dy);
  return ret;
}

std::string htos(int k) {
  std::string ret = "";
  if (k == 0) {
    ret = "u";
  } else if (k == 1) {
    ret = "d";
  } else if (k == 2) {
    ret = "l";
  } else {
    ret = "r";
  }
  return ret;
}

void* init_lowest(void* arg) { // can be threaded
  int id = *(int *)arg;
  int start = id * (N / NUM_THREADS);
  int end = (id + 1) * (N / NUM_THREADS);
  for (int i = start; i < end; i++) {
    for (int j = 0; j < N; j++) {
      int lowest = land[i][j];
      std::vector<int> heights = get_neighbors_heights(i, j);
      for (int h : heights) {
	if (h != -1 && h < lowest) {
	  lowest = h;
	}
      }
      std::string ret = "";
      if (lowest != land[i][j]) {
	for (int k = 0; k < 4; k++) {
	  if (heights[k] >= 0 && heights[k] == lowest) {
	    ret += htos(k);
	  }
	}
      }
      lowest_neighbours[i][j] = ret;
    }
  }
}

void calc_flux(int row, int col) { // can be threaded
  std::string lowest = lowest_neighbours[row][col];
  double flux_vol = 1.0;
  if (trap_water[row][col] != 0 && (lowest.length() != 0)) {
    if (trap_water[row][col] < flux_vol) {
      flux_vol = trap_water[row][col];
    }
  } else {
    flux_vol = 0.0;
  }
  flux[row][col] = flux_vol;
} 

/*void water_flow(int row, int col) {
  std::string lowest = lowest_neighbours[row][col];
  if (trap_water[row][col] > 0 && lowest.length() != 0) {
    trap_water[row][col] -= flux[row][col];
    double num = (double) lowest.length();
    for (int n = 0; n < lowest.length(); n++) {
      char direction = lowest.at(n);
      std::vector<int> dir = find_direction(direction);
      int dx = dir[0];
      int dy = dir[1];
      trap_water[row + dx][col + dy] += flux[row][col] / num;
    }
  }
  }*/
void* water_flow(void* arg) {
  int id = *(int *)arg;
  int start = id * (N / NUM_THREADS);
  int end = (id + 1) * (N / NUM_THREADS);
  for (int row = start + 1; row < end - 1; ++row) {
    for (int col = 0; col < N; ++col) {
      std::string lowest = lowest_neighbours[row][col];
      if (trap_water[row][col] > 0 && lowest.length() != 0) {
        trap_water[row][col] -= flux[row][col];
        double num = (double) lowest.length();
        for (int n = 0; n < lowest.length(); n++) {
          char direction = lowest.at(n);
	  std::vector<int> dir = find_direction(direction);
          int dx = dir[0];
          int dy = dir[1];
          trap_water[row + dx][col + dy] += flux[row][col] / num;
        }
      }
    }
  }
}

void water_flow_helper() {
  for (int thread = 0; thread < NUM_THREADS; thread++) {
    for (int head = 0; head < 2; head++) {
      for (int col = 0; col < N; col++) {
        int row;
        if (head == 0) {
          row =  thread * N / NUM_THREADS;
        }
        else {
          row = (thread + 1) * N / NUM_THREADS - 1;
        }
	std::string lowest = lowest_neighbours[row][col];
        if (trap_water[row][col] > 0 && lowest.length() != 0) {
          trap_water[row][col] -= flux[row][col];
          double num = (double) lowest.length();
          for (int n = 0; n < lowest.length(); n++) {
            char direction = lowest.at(n);
	    std::vector<int> dir = find_direction(direction);
            int dx = dir[0];
            int dy = dir[1];
            trap_water[row + dx][col + dy] += flux[row][col] / num;
          }
        }
      }
    }
  }
}

void* rain_absorb(void* arg) { // threaded                                                                                                                                          
  int id = *(int *)arg;
  int start = id * (N / NUM_THREADS);
  int end = (id + 1) * (N / NUM_THREADS);
  for (int row = start; row < end; ++row) {
    for (int col = 0; col < N; ++col) {
      if (t <= rain_time) {
	trap_water[row][col] += 1;
      }
      double water = trap_water[row][col];
      if (water > 0) {
	double remaining = water - absorb_rate;
	if (remaining < 0) {
	  absorbed[row][col] += water;
	  trap_water[row][col] = 0;
	} else {
	  absorbed[row][col] += absorb_rate;
	  trap_water[row][col] = remaining;
	}
      }
    }
  }
}

void* calculate_flux(void* arg) {
  int id = *(int *)arg;
  int start = id * (N / NUM_THREADS);
  int end = (id + 1) * (N / NUM_THREADS);
  for (int row = start; row < end; ++row) {
    for (int col = 0; col < N; ++col) {                                                                                                                                           
      calc_flux(row, col); // 这个函数可以pthread                                                                                                                                 
    }                                                                                                                                                                             
  }
}

int rain(int timestep) {
  pthread_t *threads;
  while (t <= rain_time || !dry()) {
    wet = false;
    //pthread_t *threads;
    threads = (pthread_t *) malloc(NUM_THREADS * sizeof(pthread_t));
    for (int i = 0; i < NUM_THREADS; i++) {    
      int *p = (int *) malloc(sizeof(int));  
      *p = i;    
      pthread_create(&threads[i], NULL, rain_absorb, (void *)(p));  
    }
    for (int i = 0; i < NUM_THREADS; i++) {    
      pthread_join(threads[i], NULL);  
    }
    threads = (pthread_t *) malloc(NUM_THREADS * sizeof(pthread_t));
    for (int i = 0; i < NUM_THREADS; i++) {
      int *p = (int *) malloc(sizeof(int));
      *p = i;
      pthread_create(&threads[i], NULL, calculate_flux, (void *)(p));
    }
    for (int i = 0; i < NUM_THREADS; i++) {
      pthread_join(threads[i], NULL);
    }
    /*for (int row = 0; row < N; ++row) {
      for (int col = 0; col < N; ++col) {
	water_flow(row, col); // 这个函数会有冲突
      }
      }*/

    threads = (pthread_t *) malloc(NUM_THREADS * sizeof(pthread_t));
    for (int i = 0; i < NUM_THREADS; i++) {
      int *p = (int *) malloc(sizeof(int));
      *p = i;
      pthread_create(&threads[i], NULL, water_flow, (void *)(p));
    }
    for (int i = 0; i < NUM_THREADS; i++) {
      pthread_join(threads[i], NULL);
    }
    water_flow_helper();
    if (rain_time >= t) {
      ++t;
    }
    //pthread_t *threads;
    /*threads = (pthread_t *) malloc(NUM_THREADS * sizeof(pthread_t));
    for (int i = 0; i < NUM_THREADS; i++) {    
      int *p = (int *) malloc(sizeof(int));  
      *p = i;    
      pthread_create(&threads[i], NULL, check_dry, (void *)(p));  
    }
    for (int i = 0; i < NUM_THREADS; i++) {    
      pthread_join(threads[i], NULL);  
      }*/
    ++timestep;
  }
  return timestep;
}




int main (int argc, char const *argv[]) { 
  struct timespec start_time, end_time;
  try {
    rain_time = std::stoi(argv[1]);
    absorb_rate = std::stod(argv[2]);
    N = std::stoi(argv[3]);
  }
  catch (std::exception& e) 
    {
      std::cout << e.what() << std::endl;
    }
  std::string file_name = argv[4];
  std::ifstream file(file_name);
  std::string line = "";
  init();
  int row = 0;
  while (!file.eof()) {
    std::getline(file, line);
    int col = 0;
    std::string buf;
    std::stringstream ss(line);
    while (ss >> buf) {
      try 
	{
	  land[row][col] = std::stoi(buf);
	}
      catch (std::exception e) 
	{
	  std::cout << e.what() << std::endl;
	}
      ++col;
    }
    ++row;
  }
  pthread_t *threads;
  threads = (pthread_t *) malloc(NUM_THREADS * sizeof(pthread_t));

  for (int i = 0; i < NUM_THREADS; i++) {    
    int *p = (int *) malloc(sizeof(int));  
    *p = i;    
    pthread_create(&threads[i], NULL, init_lowest, (void *)(p));  
  }
  for (int i = 0; i < NUM_THREADS; i++) {    
    pthread_join(threads[i], NULL);  
  }
  //init_lowest();
  int timestep = 0;
  clock_gettime(CLOCK_MONOTONIC, &start_time);
  int total_timestep = rain(timestep);
  clock_gettime(CLOCK_MONOTONIC, &end_time);
  double elapsed_s = calc_time(start_time, end_time) / 1000000000.0;
  std::cout << "running time is "<< elapsed_s << std::endl;
  std::cout << "Rainfall simulation took "  << total_timestep << " time steps to complete." << std::endl;
  std::cout << "The following grid shows the number of raindrops absorbed at each point: " << std::endl;
  for (std::vector<double> row : absorbed) {
    for (double value : row) {
      std::cout << std::setw(6) << std::right << value << " ";
    }
    std::cout << std::endl;
  }
  return 0;
}
