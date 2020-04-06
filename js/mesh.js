let Lx = 5; // x-width of membrane (m)
let Ly = 5; // y-width of membrane (m)
let c = 2540; // speed of sound in Mylar (membrane material) (m/s)
let dur = 0.5; // duration of simulation (seconds)
let fs = 44100; // sampling rate (Hz)
let Ns = dur * fs; // duration of simulation (samples)
let Ts = 1 / fs; // sampling period (sec)
let X = c * Ts; // distance traveled in one sampling period (m)
let Njx = Math.round(Lx / X - 1); // number of junctions along x
let Njy = Math.round(Ly / X - 1); // number of junctions along y
console.log(Njx)
console.log(Njy)
let xi = Math.ceil(Njx / 2); // junction along x receiving input excitation
let yi = Math.ceil(Njy / 2); // junction along y receiving input excitation

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Define input signal
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
let inputSignal = Array.from(Array(Ns), () => 0)
inputSignal[0] = 0.5
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Set model parameters
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
let Rw1 = 0.9; // inversion coefficient for west side of membrane
let Re2 = -0.9; // inversion coefficient for east side of membrane
let Rn3 = 0.9; // inversion coefficient for north side of membrane
let Rs4 = -0.9; // inversion coefficient for south side of membrane


let w1_end = Array.from(Array(Njy), () => 0); // value of signal at west end of membrane
let e2_end = Array.from(Array(Njy), () => 0); // value of signal at east end of membrane
let n3_end = Array.from(Array(Njy), () => 0); // value of signal at north end of membrane
let s4_end = Array.from(Array(Njy), () => 0); // value of signal at south end of membrane

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialize input and output ports. Note that the dimensions of each
// matrix are "number of junctions" by "number of rails".
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function makeNDimArray(dimensions)
{
  let t, i = 0,
    s = dimensions[0],
    arr = new Array(s);
  if (dimensions.length < 3)
  {
    for (t = dimensions[1]; i < s;)
    {
      arr[i] = new Array(t);
      arr[i].fill(0);
      i++;
    }
  }
  else
  {
    for (t = dimensions.slice(1); i < s;) arr[i++] = makeNDimArray(t);
  }
  return arr;
}

let portins = makeNDimArray([Njx, Njy, 4]);
let portouts = makeNDimArray([Njx, Njy, 4]);

for (let m = 0; m < Njx; m++)
{
  for (let q = 0; q < Njy; q++)
  {
    for (let r = 0; r < 4; r++)
    {
      portins[m][q][r] =  0;
      portouts[m][q][r] =  0;
    }
  }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Compute output signal
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
let input = 0.5;
function updateState()
{
  for (let m = 0; m < Njx; m++)
  {
    for (let q = 0; q < Njy; q++)
    {
      // Update input ports on the UPPER HORIZONTAL rails. Note that
      // the input ports on the upper rail are ahead one junction and
      // ahead one time sample compared to the output ports.
      if (m == 0)
      {
        portins[m][q][1] = w1_end[q] * Rw1;
      }
      else if (m == xi && q == yi)
      {
        portins[m][q][0] = input;
        input = 0;
      }
      else
      {
        try
        {
          portins[m][q][0] = portouts[m - 1][q][1];
        }
        catch (error)
        {
          // console.error(portouts);
          console.error(portouts[m - 1][q][1]);
        }

      }
      // Update input ports on the LOWER HORIZONTAL rails. Note that
      // the input ports on the lower rail are behind one junction and
      // ahead one time sample compared to the output ports.
      if (m == Njx - 1)
      {
        portins[m][q][1] = e2_end[q] * Re2;
      }
      else
      {
        portins[m][q][1] = portouts[m + 1][q][0];
      }


      // Update input ports on the RIGHT VERTICAL rails. Note that the
      // input ports on the right rail are ahead one junction and
      // ahead one time sample compared to the output ports.
      if (q == 0)
      {
        portins[m][q][2] = n3_end[m] * Rn3;
      }
      else
      {
        portins[m][q][2] = portouts[m][q - 1][3];
      }

      // Update input ports on the LEFT VERTICAL rails. Note that the
      // input ports on the left rail are behind one junction and
      // ahead one time sample compared to the output ports.
      if (q == Njy - 1)
      {
        portins[m][q][3] = s4_end[m] * Rs4;
      }
      else
      {
        portins[m][q][3] = portouts[m][q + 1][2];
      }

      // Set the end-point values of the membrane. Note that these
      // output port values correspond to one time sample in the past.
      w1_end[q] = portouts[1][q][0];
      e2_end[q] = portouts[Njx - 1][q][1];
      n3_end[m] = portouts[m][1][2];
      s4_end[m] = portouts[m][Njy - 1][3];
    }
  }

  for (let m = 0; m < Njx; m++)
  {
    for (let q = 0; q < Njy; q++)
    {
      // Compute the pressure Pj at the 2-port junction. Note that the
      // admittance Y = S / (rho * c), where rho is the air density
      // and c is the speed of sound in air. Since 1 / (rho * c) is
      // constant, it can be factored out of numerator and demonator
      // and canceled out. Furthermore, if we assume equal impedance
      // at the junctions, then S1=S2=S3=S4 and therefore S can be
      // factored out and cancelled.
      let Pj =
        (
          portins[m][q][0] +
          portins[m][q][1] +
          portins[m][q][2] +
          portins[m][q][3]
        ) / 2;

      // Update output ports. Note that these values will be used in
      // the next iteration of the loop because the output ports are
      // behind the input ports by one sample.
      portouts[m][q][0] = Pj - portins[m][q][0];
      portouts[m][q][1] = Pj - portins[m][q][1];
      portouts[m][q][2] = Pj - portins[m][q][2];
      portouts[m][q][3] = Pj - portins[m][q][3];
    }
  }
}
