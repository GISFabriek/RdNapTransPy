# ***********************************************************************
# Author           : Willem A. Ligtendag, De GISFabriek
# Created          : 07-14-2019
#
# Last Modified By : Willem A. Ligtendag, De GISFabriek
# Last Modified On : 07-15-2019
# ***********************************************************************
# Python PORT from C version of RDNAPTRANS
# ***********************************************************************
import transformer


def main():
    phi = 53.160753042
    lbda = 4.824761912
    h = 42.8614
    x, y, z = transformer.etrs2rdnap(phi, lbda, h)
    print("Input: Phi: " + str(phi) + "; Lambda: " + str(lbda) + "; H: " + str(h))
    print("Output: X: " + str(x) + "; Y: " + str(y) + "; Z: " + str(z))
    x = 117380.1200
    y = 575040.3400
    z = 1.0000
    phi, lbda, h = transformer.rdnap2etrs(x, y, z)
    print("Input: X: " + str(x) + "; Y: " + str(y) + "; Z: " + str(z))
    print("Output: Phi: " + str(phi) + "; Lambda: " + str(lbda) + "; H: " + str(h))


main()

