import numpy as np
import matplotlib.pyplot as plt

class Ising2D:
    
    #This function initializes the array spin_array with spins of 1 or -1.
    def __init__(self, temp=3., size=10):
        self.spin_array = np.zeros(shape=(size, size))
        for i in range(size):
            for j in range(size):
                self.spin_array[i,j] = np.random.randint(2) * 2 - 1
        #print(self.spin_array)
        
        self.temp = temp
        self.size = size
        
    #This function finds the adjacent locations for a given (x,y) location in the lattice.
    def adjacent_locs(self, x, y):
        left = [(x-1) % self.size, y]
        right = [(x+1) % self.size, y]
        down = [x, (y-1) % self.size]
        up = [x, (y+1) % self.size]
        return left, right, down, up
    
    #This function finds the energy difference between the final and potential spins in a potential spin flip.
    def energy_diff(initial_spin, final_spin, adjacent_vals):
        sum_adjacent = sum(adjacent_vals)
        energy_init = initial_spin * sum_adjacent * -1
        energy_final = final_spin * sum_adjacent * -1
        energy_diff = energy_final - energy_init
        return energy_diff
    
    
    #This function brings the system to equilibrium in the immediate cooling scenario, for flips_per_site average flips per site.    
    def equilibrium(self, flips_per_site=100):
        total_flips = flips_per_site * self.size * self.size
        for flip in range(total_flips):
            flip_x = np.random.randint(self.size)
            flip_y = np.random.randint(self.size)
            initial_spin = self.spin_array[flip_x, flip_y]
            adjacent = self.adjacent_locs(flip_x, flip_y)
            adjacent_vals = [self.spin_array[adjacent_loc[0],adjacent_loc[1]] for adjacent_loc in adjacent]
            delta_E = Ising2D.energy_diff(initial_spin, -1 * initial_spin, adjacent_vals)
            if delta_E <= 0:
                self.spin_array[flip_x, flip_y] *= -1
            else:
                rand = np.random.random_sample()
                prob = np.exp(-1 * delta_E / self.temp)
                if rand < prob:
                    self.spin_array[flip_x, flip_y] *= -1
    
    
    #This function simulates cooling the system gradually to low temperature with a characteristic time of cooling_time for 
    #flips_per_site average flips per site.
    def cool(self, temp_start, flips_per_site=1000, cooling_time=100.):
        for n in range(flips_per_site):
            self.temp = temp_start * np.exp(-1 * n / cooling_time)
            for i in range(self.size * self.size):
                flip_x = np.random.randint(self.size)
                flip_y = np.random.randint(self.size)
                initial_spin = self.spin_array[flip_x, flip_y]
                adjacent = self.adjacent_locs(flip_x, flip_y)
                adjacent_vals = [self.spin_array[adjacent_loc[0],adjacent_loc[1]] for adjacent_loc in adjacent]
                delta_E = Ising2D.energy_diff(initial_spin, -1 * initial_spin, adjacent_vals)
                if delta_E <= 0:
                    self.spin_array[flip_x, flip_y] *= -1
                else:
                    rand = np.random.random_sample()
                    prob = np.exp(-1 * delta_E / self.temp)
                    if rand < prob:
                        self.spin_array[flip_x, flip_y] *= -1
                
            
                
    #This function finds the average energy per site of the lattice.        
    def get_average_energy_per_site(self):
        double_energy = 0
        for i in range(self.size):
            for j in range(self.size):
                spin_ij = self.spin_array[i,j]
                adjacent = self.adjacent_locs(i, j)
                adjacent_vals = [self.spin_array[adjacent_loc[0],adjacent_loc[1]] for adjacent_loc in adjacent]
                sum_adjacent = sum(adjacent_vals)
                energy_ij = spin_ij * sum_adjacent * -1
                double_energy += energy_ij
        return double_energy / (2.0 * self.size * self.size)
    
    #This function finds the average spin of the states in the lattice
    def get_average_spin_per_site(self):
        total_spin = 0
        for i in range(self.size):
            for j in range(self.size):
                spin_ij = self.spin_array[i,j]
                total_spin += spin_ij
        return total_spin / (self.size * self.size)
    
    #This function displays the lattice pattern with yellow as spin up and blue as spin down.
    def display(self, title = ''):
        plt.figure()
        plt.imshow(self.spin_array, vmin = -1, vmax = 1)
        plt.title(title)
        
    #This function returns the string representation of the temperature, average energy, and average spin of the system.
    def __str__(self):
        return f'Temperature: {self.temp}, Average Energy: {self.get_average_energy_per_site()}, Average Spin: {self.get_average_spin_per_site()}'
    
    
lattice = Ising2D(temp=3., size=10)
lattice.display(title='Start')

#These graphs are the equilibrium states of a two-state paramagnet 
# that is immediately cooled to zero temperature
for i in range(12):
    lattice = Ising2D(temp=3., size=10) # Initialize 
    lattice.temp = 0.001      # Rapid quenching to low T
    lattice.equilibrium(100) # Allow to come to equilibrium
    print(lattice)            # Print temp, energy and spin per site
    
    # Display graphically
    title = 'Rapid T=0.001 Spin:{:5.2f}'.format(lattice.get_average_spin_per_site())
    title += ' Energy:{:5.2f}'.format(lattice.get_average_energy_per_site())
    lattice.display(title=title)

#These graphs are the equilibrium states of a two-state paramagnet 
# that is gradually cooled to zero temperature

for i in range(10):
    lattice = Ising2D(temp=3., size=10)
    lattice.cool(temp_start=3., flips_per_site=500, cooling_time=100)
    print(lattice)
    
    # Display graphically
    title = 'Slow T=0.02 Spin:{:5.2f}'.format(lattice.get_average_spin_per_site())
    title += ' Energy:{:5.2f}'.format(lattice.get_average_energy_per_site())
    lattice.display(title=title)