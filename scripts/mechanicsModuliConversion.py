import argparse

def computeModuli_from_Young_poisson(E, nu):
  K = E / (3.0 * (1.0 - 2.0 * nu))
  G = E / (2.0 * (1.0 + nu))
  return K, G

def computeModuli_from_bulk_shear(K, G):
   nu = ( 3.0 * K - 2.0 * G ) / (2.0*G + 6.0 * K)
   E  = G * ( 2.0 * (1.0 + nu) )
   return E, nu

def computeModuli_from_bulk_poisson(K, nu):
   E = 3 * K * ( 1 - 2 * nu) 
   G  = E / (2.0 * (1.0 + nu))
   return E, G

def computeModuli_from_shear_poisson(G, nu):
   E = 2 * G * (1+nu)
   K =  E / (3.0 * (1.0 - 2.0 * nu))
   return E, K   

def parseArguments():
    """
    Parse command line arguments into an ArgumentParser instance.
    :param arguments: The command line arguments as an array of string.
    :return: Values, inputType.
    """
    parser = argparse.ArgumentParser(description=".")

    parser.add_argument('m1', 
                         metavar='M1', 
                         type=float,
                         help='The first mechanical modulus')

    parser.add_argument('m2', 
                         metavar='M2', 
                         type=float,
                         help='The second mechanical modulus')                     

    parser.add_argument("-i",
                        "--input",
                        dest="input_type", 
                        type=str,
                        choices=["YoungPoisson", "BulkShear", "BulkPoisson", "ShearPoisson"],
                        default="YoungPoisson",
                        help="Name of the pair of mechanical moduli provided.")                 

    args = parser.parse_args()

    Values = [args.m1, args.m2]

    inputModuliType = str(args.input_type)

    return Values, inputModuliType

def convert(values, inputModuliType):
    
    if (inputModuliType == "YoungPoisson"):
        E = values[0]
        nu = values[1]
        K, G = computeModuli_from_Young_poisson(E, nu)
    elif (inputModuliType == "BulkShear"): 
        K = values[0]
        G = values[1]
        E, nu = computeModuli_from_bulk_shear(K, G)
    elif (inputModuliType == "BulkPoisson"):
        K = values[0]
        nu = values[1]
        E, G = computeModuli_from_bulk_poisson(K, nu)
    elif (inputModuliType == "ShearPoisson"):
        G = values[0]
        nu = values[1]
        E, K = computeModuli_from_shear_poisson(K, G)

    print("Mechanics Moduli:")
    print("Young's modulus = {} ".format(E))
    print("Bulk modulus = {} ".format(K))
    print("Shear modulus = {} ".format(G))
    print("Poisson ratio = {} ".format(nu))
    

if __name__ == '__main__':    
    Values, inputModuliType = parseArguments() 
    convert(Values, inputModuliType)
