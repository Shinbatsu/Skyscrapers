#include <vector>
#include <stack>
#include <iostream>



#define row_index(x, N) ((x > (N - 1) && x < N * 2) * (x - N) + (x > (N * 2 - 1) && x < N * 3) * (N - 1) + (x > (N * 3 - 1)) * (N * 4 - 1 - x))
#define col_index(x, N) ((x < N) * (x) + (x > (N - 1) && x < (2 * N)) * (N - 1) + (x > (2 * N - 1) && x < (3 * N)) * ((3 * N - 1) - x))


std::vector<std::vector<int>> copy_vector(std::vector<std::vector<int>> scrapers, int N) {
  std::vector<std::vector<int>> vec(N);
  for(int i = 0; i < N; i++) {
    for(auto plot: scrapers.at(i)) vec.at(i).push_back(plot);
  }
  return vec;
}



int degrees_of_freedom(int plot, int N) {
  int dof = 0;
  for(int i = 0; i < N; i++) {
    dof += (plot >> i) & 1;
  }
  return dof;
}


void generate_vectors(std::vector<std::vector<int>> &scrapers, std::vector<std::vector<int*>> &vectors, int i, int N) {
  for(int x = 0; x < N; x++) {
    vectors.at(0).at(x) = (&(*(scrapers.at(x).begin() + i)));
    vectors.at(1).at(x) = (&(*(scrapers.at(i).begin() + x)));
    vectors.at(2).at(x) = (&(*(scrapers.at(N - 1 - x).begin() + i)));
    vectors.at(3).at(x) = (&(*(scrapers.at(i).begin() + N - 1 - x)));
  }
}
void placeScraper(std::vector<std::vector<int>> &scrapers, int row, int col, int target, int N) {
  scrapers.at(row).at(col) = 0;
  for(int purge_index = 0; purge_index < N; purge_index++) {
    
    int dof = degrees_of_freedom(scrapers.at(row).at(purge_index), N);
    scrapers.at(row).at(purge_index) &= ~target;
    if(dof == 2 && degrees_of_freedom(scrapers.at(row).at(purge_index), N) == 1) {
      placeScraper(scrapers, row, purge_index, scrapers.at(row).at(purge_index), N);
    }
    
    dof = degrees_of_freedom(scrapers.at(purge_index).at(col), N);
    scrapers.at(purge_index).at(col) &= ~target;
    if(dof == 2 && degrees_of_freedom(scrapers.at(purge_index).at(col), N) == 1) {
      placeScraper(scrapers, purge_index, col, scrapers.at(purge_index).at(col), N);
    }
  }
  
  scrapers.at(row).at(col) = target;
}
void preprocess(std::vector<std::vector<int>> &scrapers, std::vector<std::vector<int>> starting_grid, const std::vector<int> &clues, int N) {

  std::vector<std::vector<int *>> vectors(4, std::vector<int *>(N));
  int i;

  for(i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      if(starting_grid.at(i).at(j) != 0) placeScraper(scrapers, i, j, (1 << (starting_grid.at(i).at(j) - 1)), N);
    }
  }
  
  for(i = 0; i < N; i++) {
    int bearings[] = {clues[i], clues[N * 4 - 1 - i], clues[N * 3 - 1 - i], clues[i + N]};
    generate_vectors(scrapers, vectors, i, N);
    
    for(int x = 0; x < 4; x++) {
        
      int check_list = 0;
      for(int j = 0; j < bearings[x] - 1; j++) {
          
        check_list += 1 << (N - 1 - j);
        *vectors.at(x).at(bearings[x] - 2 - j) &= ~check_list;
        
      }
    }
  }
  
  for(i = 0; i < N * 4; i++) {
    if(clues.at(i) == 1) placeScraper(scrapers, row_index(i, N), col_index(i, N), 1 << (N - 1), N);
  }
}
std::vector<std::vector<int>> get_output_vector(std::vector<std::vector<int>> scrapers, int N) {
  std::vector<std::vector<int>> full_output(N);
  for(int i = 0; i < N; i++) {
    for(int k = 0; k < N; k++) {
      for(int j = 0; j < N; j++) {
        scrapers.at(i).at(k) += (j + 1 - scrapers.at(i).at(k)) * (scrapers.at(i).at(k) == 1 << j);
      }
      full_output.at(i).push_back(scrapers.at(i).at(k));
    }
  }
  return full_output;
}


void resolve_only_possibilities(std::vector<std::vector<int>> &scrapers, std::vector<std::vector<int*>> &vectors, int i, int N) {
    for(int x = 0; x < 2; x++) {
      for(int target = 0; target < N; target++) {
        std::vector<int> check;
        for(int location = 0; location < N; location++) {
          if(*vectors.at(x).at(location) & 1 << target) check.push_back(location);
        }
        
        if(check.size() == 1 && (*vectors.at(x).at(check.back()) & ~(1 << target))) {
          int col = x * ((N - 1 - check.back()) * (x >> 1) + check.back() * !(x >> 1)) + !x * i;
          int row = !x * ((N - 1 - check.back()) * (x >> 1) + check.back() * !(x >> 1)) + x * i;
          placeScraper(scrapers, row, col, (1 << target), N);
        }
      }
    }
}
std::vector<std::vector<int>> SolvePuzzle(const std::vector<int> &clues,std::vector<std::vector<int>> starting_grid, int N) {
  std::vector<std::vector<std::vector<int>>> solution_vector;
  std::stack<std::vector<std::vector<int>>> scraper_stack;
  std::vector<std::vector<int>> scrapers(N, std::vector<int>(N, (1 << N) - 1));
  std::vector<std::vector<int *>> vectors(4, std::vector<int *>(N));
  int i;
  
  preprocess(scrapers, starting_grid, clues, N);
  int last_dof = 0;
  
  while(solution_vector.size() < 1) {
    for(i = 0; i < N; i++) {
        generate_vectors(scrapers, vectors, i, N);
        resolve_only_possibilities(scrapers, vectors, i, N);
    }
    
    
    for(i = 0; i < N; i++) {
        
      generate_vectors(scrapers, vectors, i, N);
      int bearings[] = {clues[i], clues[N * 4 - 1 - i], clues[N * 3 - 1 - i], clues[i + N]};
     
      for(int x = 0; x < 2; x++) {
        int line_dof = 0;
        
        for(int k = 0; k < N; k++) {
          line_dof += degrees_of_freedom(*vectors.at(x).at(k), N);
        }
        if(line_dof > (N * N) || line_dof == N || (bearings[x] == 0 && bearings[x + 2] == 0)) continue;
      
        std::stack<std::vector<int>> p_space;
        std::vector<std::vector<int>> completed;
        std::vector<int> solution_vector;
      
        for(int k = 0; k < N; k++) {
          solution_vector.push_back(*vectors.at(x).at(k));
        }
      
        p_space.push(solution_vector);
        while(p_space.size()) {
          std::vector<int> line = p_space.top();
          p_space.pop();
        
          int dof, k = 0;
          for(k = 0; k < N; k++) {
            dof = degrees_of_freedom(line.at(k), N);
            if(dof != 1) break;
          }
        
          if(dof == 0) continue;
        
          if(k == N) {
            int watermarkLeft = line.at(0);
            int watermarkRight = line.at(N - 1);
            int risesLeft = 1;
            int risesRight = 1;
            for(int j = 1; j < N; j++) {
              if(line.at(j) > watermarkLeft) {
                watermarkLeft = line.at(j);
                risesLeft++;
              }
              if(line.at(N - 1 - j) > watermarkRight) {
                watermarkRight = line.at(N - 1 - j);
                risesRight++;
              }
            }
            if((risesLeft == bearings[x] || bearings[x] == 0) && (risesRight == bearings[x + 2] || bearings[x + 2] == 0)) completed.push_back(line);
            continue;
          }
          for(int j = 0; j < N; j++) {
            if(~line.at(k) & (1 << j)) continue;
            
            std::vector<int> temp;
            for(int m = 0; m < N; m++) {
              temp.push_back(line.at(m) & ~(1 << j));
            }
            temp.at(k) = 1 << j;
            p_space.push(temp);
          }
        }
        if(completed.size()) {
          while(completed.size() > 1) {
            for(int k = 0; k < N; k++) {
              completed.at(0).at(k) |= completed.back().at(k);
            }
            completed.pop_back();
          }


          for(int k = 0; k < N; k++) {
            if(degrees_of_freedom(completed.back().at(k), N) == 1 && completed.back().at(k) != *vectors[x][k]) {
              int col = x * ((N - 1 - k) * (x >> 1) + k * !(x >> 1)) + !x * i;
              int row = !x * ((N - 1 - k) * (x >> 1) + k * !(x >> 1)) + x * i;
              placeScraper(scrapers, row, col, completed.back().at(k), N);
            } else {
              *vectors.at(x).at(k) = completed.back().at(k);
            }
            
          }
          
        } else {
          if(scraper_stack.size() == 0) break;
          last_dof = 0;
          scrapers = scraper_stack.top();
          scraper_stack.pop();
        }
      }
    }
      int total_dof = 0;
      for(int k = 0; k < N; k++) {
        for(int j = 0; j < N; j++) {
          total_dof += degrees_of_freedom(scrapers.at(k).at(j), N);
        }
      }


      if(total_dof == N * N) {
          
        int bad_vals = 0;
        for (int j = 0; j < N; j++) {
          int bearings[] = {clues.at(j), clues.at(N * 4 - 1 - j), clues.at(N * 3 - 1 - j), clues.at(j + N)};
          generate_vectors(scrapers, vectors, j, N);
          for (int x = 0; x < 2; x++) {
            int watermarkLeft = *vectors.at(x).at(0);
            int watermarkRight = *vectors.at(x).at(N - 1);
            int risesLeft = 1;
            int risesRight = 1;
            
            for(int k = 1; k < N; k++) {
              if(*vectors.at(x).at(k) > watermarkLeft) {
                watermarkLeft = *vectors[x][k];
                risesLeft++;
              }
            }
            for(int k = 1; k < N; k++) {
              if(*vectors.at(x).at(N - 1 - k) > watermarkRight) {
                watermarkRight = *vectors.at(x).at(N - 1 - k);
                risesRight++;
              }
            }
            if((risesLeft != bearings[x] && bearings[x] != 0) || (risesRight != bearings[x + 2] && bearings[x + 2] != 0)) {
              bad_vals++;
              
            }
          }
        }
        if(bad_vals == 0) {
            solution_vector.push_back(get_output_vector(scrapers, N));
        }
        if(scraper_stack.size() == 0) break;
        last_dof = 0;
        scrapers = scraper_stack.top();
        scraper_stack.pop();
      }
      if(total_dof == last_dof) {
        int row = -1, col = -1;
        for(int row_ind = 0; row_ind < N && col == -1; row_ind++) {
          for(int col_ind = 0; col_ind < N && col == -1; col_ind++) {
            if(degrees_of_freedom(scrapers.at(row_ind).at(col_ind), N) > 1) {
              row = row_ind;
              col = col_ind;
            }
          }
        }
        for(int target = 0; target < N; target++) {
          if(!(scrapers.at(row).at(col) >> target & 1)) continue;
          
          std::vector<std::vector<int>> temp = copy_vector(scrapers, N);
          placeScraper(temp, row, col, 1 << target, N);
          scraper_stack.push(temp);
        }
        
        if(scraper_stack.size() == 0) break;
        
        scrapers = scraper_stack.top();
        scraper_stack.pop();
      }
      
      last_dof = total_dof;
      
  }
  
  return solution_vector.at(0);
}