import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

def define_beam(mean, std, num_of_particles):
    return np.random.normal(mean, std, (num_of_particles, 4)).T

def F(f):
    return np.array([[1,    0, 0,   0], 
                     [-1/f, 1, 0,   0],
                     [0,    0, 1,   0],
                     [0,    0, 1/f, 1]])
def D(s):
    return np.array([[1, s, 0, 0], 
                     [0, 1, 0, 0],
                     [0, 0, 1, s],
                     [0, 0, 0, 1]]) 

def beam_line(beam, beamline_lattice):

    beam_lattice_result = []
    beam_lattice_result.append(beam)

    for magnet_matrix in beamline_lattice[-1::-1]:

        beam = magnet_matrix@beam
        beam_lattice_result.append(beam)

    return beam_lattice_result

def plot(beam_matrices):
    
    fig, axes = plt.subplots(len(beam_matrices), 2, sharex = True, figsize = (10, 14))
    counter = 0

    for beam in beam_matrices:

        x_coord = beam[0]
        x_vel   = beam[1]
        y_coord = beam[2]
        y_vel   = beam[3]

        axes[counter, 0].scatter(x = x_coord , y = x_vel ,color ="b")
        axes[counter, 0].set(xlabel = "x [mm]", ylabel = "x' [mrad]")

        axes[counter, 1].scatter(x = y_coord , y = y_vel ,color ="r")
        axes[counter, 1].set(xlabel = "y [mm]", ylabel = "y' [mrad]")

        counter += 1

    plt.show()

def main():

    """
    USER INPUTS
    """
    drift1 = float(input("Enter drift length before quadrupole: "))
    drift2 = float(input("Enter drift length after quadrupole: "))
    f      = float(input("Focal lenght of the quadrupole: "))

    beam = define_beam(0, 1, 10000)                                        # Create a 2D beam with mean = 0, standard deviation = 1, number of particles = 10,000 
    beam_lattice = [D(drift1), F(f), D(drift2), F(-f), D(10)]              # Design a beam lattice with user inputs
    beam_matrices = beam_line(beam, beam_lattice)                          # Send the particles to the beam lattice and return results
    plot(beam_matrices)                                                    # Plot the results

if __name__ == "__main__":

    main()