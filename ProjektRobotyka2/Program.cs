using System;
using System.Collections.Generic;
using System.Globalization;


namespace ProjektRobotyka2
{
    class Program
    {
        const double RadianCycle = 2 * Math.PI;

        private static double GetInputs(string msg, bool allowNegative = true)
        {
            Console.WriteLine(msg);
            double output;
            while (true)
            {
                try
                {
                    output = Double.Parse(Console.ReadLine(), CultureInfo.InvariantCulture);
                    if (!allowNegative && output <= 0) throw new FormatException();
                        break;
                }
                catch (FormatException)
                {
                    Console.WriteLine("Incorrect input format!");
                }
            }
            return output;
        }


        private static string ShowValue(double d)
        {
            return d.ToString(CultureInfo.InvariantCulture);
        }

        private static double[] CalculatePositionOfFirstJoint(double startingPos,double X,double Y, double Z)
        {    
            startingPos = ClampRadians(startingPos);
            double[] result;

            if (X == 0 && Y == 0)
            {
                result = new double[1] { startingPos };
                return result;
            }

            double option1 = Math.Atan2(Y, X);
            double option2 = option1 + Math.PI;

            result = new double[2] { option1, option2 };
            return result;
        }

        private static double[] CalculatePositionOfThirdJoint(double X, double Y, double Z, double finalTheta1)
        {

            double sub1 = Math.Pow(Math.Cos(finalTheta1) * X + Math.Sin(finalTheta1) * Y - 16, 2);

            double sub2 = Math.Pow(Z - 350, 2);

            double sub3 = 2 * 220 * 220;

            double sub = ((sub1 + sub2) / sub3) - 1;

            if (sub > 1 || sub < -1)
            { 
                return null;
            }

            double option1 = Math.Atan2(Math.Sqrt(1-(sub*sub)),sub);
            double option2 = Math.Atan2(-Math.Sqrt(1 - (sub * sub)), sub);

            return new double[2] { option1, option2 };
        }


        private static double CalculatePositionOfSecondJoint(double X, double Y, double Z, double finalTheta1,double finalTheta3)
        {
            double sub1 = 220 * (1 + Math.Cos(finalTheta3)) * (Z - 350) - (220 * Math.Sin(finalTheta3)) * (Math.Cos(finalTheta1) * X + Math.Sin(finalTheta1) * Y - 16);
            double sub2 = 220 * (1 + Math.Cos(finalTheta3)) * (Math.Cos(finalTheta1) * X + Math.Sin(finalTheta1) * Y - 16) - (220 * Math.Sin(finalTheta3)) * (Z - 350);

            double option = Math.Atan2(sub1, sub2);

            return ClampRadians(option);
        }

        private static double[] SelectCorrectOption(double startTheta1, double startTheta2, double startTheta3,List<double[]> configurations)
        {
            startTheta1 = ClampRadians(startTheta1);
            startTheta2 = ClampRadians(startTheta2);
            startTheta3 = ClampRadians(startTheta3);

            double totalDelta = 3 * Math.PI;
            double[] result= new double[3];

            for (int i = 0; i < configurations.Count; i++)
            {
                double[] d = configurations[i];
                double currentDelta = radDifference(startTheta1, d[0]);
                currentDelta += radDifference(startTheta2, d[1]);
                currentDelta += radDifference(startTheta3, d[2]);

                if (currentDelta < totalDelta)
                {
                    totalDelta = currentDelta;
                    result = d;
                }
            }
            return result;
        } 

        private static double radDifference(double r1, double r2)
        {
            double rDiff = Math.Abs(r1 - r2);
            if (rDiff > Math.PI)
            {
                rDiff = 2 * Math.PI - rDiff;
            }
            return rDiff;
        }

        private static double ClampRadians(double rad)
        {
            int tmp = (int)(rad/(RadianCycle));

            return rad-(tmp*RadianCycle);
        }



        private static double CalculateAngularAcceleration(double theta0, double theta1, double time)
        {
            double deltaTheta = radDifference(theta0 , theta1);
            double alpha = 4 * deltaTheta / (time * time);

            return alpha;
        }

        private static double CalculateMaxAngularVelocity(double alpha, double t)
        {
            double omega = (alpha * t) / 2;
            return omega;
        }

        static void Main(string[] args)
        {
            bool lock1 = true;
            while (lock1)
            {
                Console.Clear();
                char[] input = new char[0];
                Console.WriteLine("Please select which excercises do you wish to run.");
                Console.WriteLine("Write 1 for exercise 1.1 (Wyznacz model matematyczny trajektorii bang-bang)");
                Console.WriteLine("Write 2 for exercise 1.2 („Realizacja” trajektorii w rzeczywistym manipulatorze)");
                Console.WriteLine("Write 3 for exercise 2.1 (Wyznacz model matematyczny trajektorii LSPB)");
                Console.WriteLine("Write 4 for exercise 2.2 („Realizacja” trajektorii liniowej w rzeczywistym manipulatorze.)");
                Console.WriteLine("Please note that you can write multiple numbers to run multiple exercises at the same time (for example write 1234 to run all four)");
                bool correct = false;
                while (!correct)
                {
                    correct = true;
                    input = Console.ReadLine().ToCharArray();
                    foreach (char c in input)
                    {
                        if (c != '1' && c != '2' && c != '3' && c != '4')
                        {
                            correct = false;
                        }
                    }
                    if (!correct)
                    {
                        Console.WriteLine("Incorrect input. Please try again, remember to not put any separator between the numbers!");
                    }
                }
                bool[] exercisesToRun = new bool[4] { false, false, false, false };

                foreach (char c in input)
                {
                    switch (c)
                    {
                        case '1':
                            exercisesToRun[0] = true;
                            break;
                        case '2':
                            exercisesToRun[1] = true;
                            break;
                        case '3':
                            exercisesToRun[2] = true;
                            break;
                        case '4':
                            exercisesToRun[3] = true;
                            break;
                    }
                }


                if (exercisesToRun[0])
                {
                    if (!Part1Main(exercisesToRun[1]))
                    {
                        Console.WriteLine("The chosen trajectory is impossible to achive by this robot.");
                    }
                }


                if (exercisesToRun[2])
                {
                    if (!Part2Main(exercisesToRun[3]))
                    {
                        Console.WriteLine("The chosen trajectory is impossible to achive by this robot.");
                    }
                }
                Console.WriteLine();
                bool lock2 = true;
                while (lock2)
                {
                    Console.WriteLine("Do you want to restart the application (y/n)?");
                    char key = Console.ReadKey().KeyChar;
                    if (key=='y')
                    {
                        lock2 = false;
                    }
                    else if (key== 'n')
                    {
                        lock2 = false;
                        lock1 = false;
                    }
                }
             }
            Console.ReadKey();
        }

        static bool Part1Main(bool runSecondPart)
        {
            Console.WriteLine();
            Console.Write("Data for exercise 1.1");
            if (runSecondPart) Console.WriteLine(" and 1.2");
            Console.WriteLine();
            Console.WriteLine("Assume that starting velocity is zero (v0=0) and starting time is zero (t0=0).");
            Console.WriteLine("Please remember that all input values are treated as real numbers. Remember to use the dot (.) as a decimal separator!");
            Console.WriteLine();

            double theta1 = GetInputs("Input the starting position of the first joint in radians (theta1): ");
            double theta2 = GetInputs("Input the starting position of the second joint in radians (theta2): ");
            double theta3 = GetInputs("Input the starting position of the third joint in radians (theta3): ");

            double positionX = GetInputs("Input the target x coordinate of the robot's manipulator (in mm): ");
            double positionY = GetInputs("Input the target y coordinate of the robot's manipulator (in mm): ");
            double positionZ = GetInputs("Input the target z coordinate of the robot's manipulator (in mm): ");

            double t1 = GetInputs("Input the target time of the trajectory, provided that the starting time is t0 = 0 seconds (in seconds): ",false);

            List<double[]> finalConfigurations = Part1CalcConfigs(theta1, positionX, positionY, positionZ);
            if (finalConfigurations.Count == 0)
            {
                return false;
            }

            double[] optimalFinalConfig = SelectCorrectOption(theta1, theta2, theta3, finalConfigurations);

            

            double alpha1 = CalculateAngularAcceleration(theta1, optimalFinalConfig[0], t1);
            double alpha2 = CalculateAngularAcceleration(theta2, optimalFinalConfig[1], t1);
            double alpha3 = CalculateAngularAcceleration(theta3, optimalFinalConfig[2], t1);

            double omega1 = CalculateMaxAngularVelocity(alpha1, t1);
            double omega2 = CalculateMaxAngularVelocity(alpha2, t1);
            double omega3 = CalculateMaxAngularVelocity(alpha3, t1);

            Console.WriteLine();
            Console.WriteLine("Results (1.1): ");
            Console.WriteLine();
            Console.WriteLine("The target joint positions ([theta1(t1),theta2(t1),theta3(t3)]) = ["+ShowValue(optimalFinalConfig[0])+", "+ ShowValue(optimalFinalConfig[1]) + ", " + ShowValue(optimalFinalConfig[2])+"]");
            Console.WriteLine("The angular acceleration values ([theta1(t1),theta2(t1),theta3(t3)]) = [" + ShowValue(alpha1) + ", " + ShowValue(alpha2) + ", " + ShowValue(alpha3) + "]");
            Console.WriteLine("The highest achivable angular velocities ([theta1(t1),theta2(t1),theta3(t3)]) = [" + ShowValue(omega1) + ", " + ShowValue(omega2) + ", " + ShowValue(omega3) + "]");

            if (runSecondPart)
            {
                Console.WriteLine();
                double deltaT = GetInputs("Input the value of delta t, for the discreet representation of the trajectory (1.2) (in seconds): ", false);

                double[] theta1OverTime = Part1Point2(optimalFinalConfig[0], t1, alpha1, deltaT,omega1);
                double[] theta2OverTime = Part1Point2(optimalFinalConfig[1], t1, alpha2, deltaT,omega2);
                double[] theta3OverTime = Part1Point2(optimalFinalConfig[2], t1, alpha3, deltaT,omega3);

                Console.WriteLine();
                Console.WriteLine("Results (1.2): ");
                Console.WriteLine();

                Console.Write("The list of consecutive positions for joint 1 ([theta1(t0), theta1(t0+deltat)+...+theta1(t1)]): [" + ShowValue(theta1OverTime[0]));
                for (int i = 1; i < theta1OverTime.Length; i++)
                {
                    Console.Write(", " +ShowValue(theta1OverTime[i]));
                }
                Console.WriteLine("]");

                Console.Write("The list of consecutive positions for joint 2 ([theta2(t0), theta2(t0+deltat)+...+theta2(t1)]): [" + ShowValue(theta2OverTime[0]));
                for (int i = 1; i < theta2OverTime.Length; i++)
                {
                    Console.Write(", " + ShowValue(theta2OverTime[i]));
                }
                Console.WriteLine("]");

                Console.Write("The list of consecutive positions for joint 3 ([theta3(t0), theta3(t0+deltat)+...+theta3(t1)]): [" + ShowValue(theta3OverTime[0]));
                for (int i = 1; i < theta3OverTime.Length; i++)
                {
                    Console.Write(", " + ShowValue(theta3OverTime[i]));
                }
                Console.WriteLine("]");
            }



            return true;
        }
        static List<double[]> Part1CalcConfigs(double theta1,double positionX, double positionY, double positionZ)
        {
            List<double[]> finalConfigurations = new List<double[]>();
            double[] finalThetas1 = CalculatePositionOfFirstJoint(theta1, positionX, positionY, positionZ);
            foreach (double d in finalThetas1)
            {

                finalConfigurations.Add(new double[3] { d, 0, 0 });
            }
            int originalCount = finalConfigurations.Count;
            for (int i = 0; i < originalCount; i++)
            {
                double[] d = finalConfigurations[i];
                double[] finalThetas3 = CalculatePositionOfThirdJoint(positionX, positionY, positionZ, d[0]);


                if (finalThetas3 == null)
                {
                    finalConfigurations.RemoveAt(i);
                    i--;
                    originalCount--;
                }
                else
                {
                    d[2] = finalThetas3[0];
                    finalConfigurations.Add(new double[3] { d[0], 0, finalThetas3[1] });
                }
            }


            for (int i = 0; i < finalConfigurations.Count; i++)
            {
                double[] d = finalConfigurations[i];
                d[1] = CalculatePositionOfSecondJoint(positionX, positionY, positionZ, d[0], d[2]);
            }
            return finalConfigurations;
        }


        public static double[] Part1Point2(double theta0, double t1, double alpha, double deltat, double Vmax)
        {

            int k = 1;
            double[] result = new double[(int)Math.Ceiling(t1 / deltat) + 1];
            result[0] = theta0+(alpha * deltat * deltat) * 0.5;
            int kSwitch = 0;
            while (k * deltat <= t1)
            {
                if (k * deltat <= t1 * 0.5)
                {
                    result[k] = result[k - 1] + (k * deltat * alpha) + (alpha * deltat * deltat * 0.5);
                }
                if (k * deltat > t1 * 0.5)
                {
                    if (kSwitch == 0) kSwitch = k - 1;
                    result[k] = result[k - 1] + (Vmax - (k - kSwitch) * deltat * alpha) - (alpha * deltat * deltat * 0.5);
                }
                k++;
            }
            return result;
        }



        static double CalculateXCoordinate(double t1, double t2, double t3)
        {
            double X = Math.Cos(t1) * (16 + 220 * ((-Math.Sin(t2) * Math.Sin(t3)) + (Math.Cos(t2) * (1 + Math.Cos(t3)))));

            return X;
        }

        static double CalculateYCoordinate(double t1, double t2, double t3)
        {
            double Y = Math.Sin(t1) * (16 + 220 * ((-Math.Sin(t2) * Math.Sin(t3)) + (Math.Cos(t2) * (1 + Math.Cos(t3)))));

            return Y;
        }

        static double CalculateZCoordinate(double t2, double t3)
        {
            double Z = 220 *(Math.Sin(t2) * (1 + Math.Cos(t3))+ Math.Sin(t3)*Math.Cos(t2))+350;

            return Z;
        }

        static double[] CalculateDirectionVector(double[] pos0, double[] pos1, double coef)
        {
            double[] k = new double[3];

            k[0] = (pos1[0] - pos0[0]) * 1 / coef;
            k[1] = (pos1[1] - pos0[1]) * 1 / coef;
            k[2] = (pos1[2] - pos0[2]) * 1 / coef;

            return k;
        }
        static double CalcCoef(double[] pos0, double[] pos1)
        {
            double coef = Math.Sqrt((pos1[0] - pos0[0])*(pos1[0] - pos0[0]) + (pos1[1] - pos0[1])*(pos1[1] - pos0[1]) + (pos1[2] - pos0[2])*(pos1[2] - pos0[2]));
            return coef;
        }

        static double CalcFirstSegmentTime(double coef,double V, double t1)
        {
            double tb = (-coef + V * t1) / V;
            return tb;
        }

        /// <summary>
        /// Calculates acceleration in mm/s^2.
        /// </summary>
        /// <param name="time">Time in seconds</param>
        /// <param name="velocity">Velocity in mm/s</param>
        /// <returns></returns>
        static double CalculateLinearAcceleration(double time, double velocity)
        {
            return velocity / time;
        }

        static double[] CalculateParametersOfFirstEquation(double acc)
        {
            double[] result = new double[3] { acc / 2, 0, 0 };
            return result;
        }
        //UNFINISHED
        static double[] CalculateParametersOfSecondEquation(double V,double coef,double t1)
        {
            double[] result = new double[3] { 0, V, (coef-V*t1)/2 };
            return result;
        }

        static double[] CalculateParametersOfThirdEquation(double coef,double acc,double t1)
        {
            double[] result = new double[3] { acc/2, acc*t1, coef - (acc / 2) * t1 * t1 };
            return result;
        }

        static bool Part2Main(bool runSecondPart)
        {
            Console.WriteLine();
            Console.Write("Data for exercise 2.1");
            if (runSecondPart) Console.WriteLine(" and 2.2");

            Console.WriteLine("Assume that starting velocity is zero (v0=0) and starting time is zero (t0=0).");
            Console.WriteLine("Please remember that all input values are treated as real numbers. Remember to use the dot (.) as a decimal separator!");
            Console.WriteLine();

            double theta1 = GetInputs("Input the starting position of the first joint in radians (theta1): ");
            double theta2 = GetInputs("Input the starting position of the second joint in radians (theta2): ");
            double theta3 = GetInputs("Input the starting position of the third joint in radians (theta3): ");

            double positionX = GetInputs("Input the target x coordinate of the robot's manipulator (in mm): ");
            double positionY = GetInputs("Input the target y coordinate of the robot's manipulator (in mm): ");
            double positionZ = GetInputs("Input the target z coordinate of the robot's manipulator (in mm): ");

            if (Part1CalcConfigs(theta1, positionX, positionY, positionZ).Count == 0)
            {
                return false;
            }


            double t1 = GetInputs("Input the target time of the trajectory, provided that the starting time is t0 = 0 seconds (in seconds): ", false);
            double V = GetInputs("Input the linear velocity V during the second segment of the trajectory (in mm/s): ", false);

            double[] position0 = new double[3] { CalculateXCoordinate(theta1, theta2, theta3), CalculateYCoordinate(theta1, theta2, theta3), CalculateZCoordinate(theta2, theta3) };
            double[] position1 = new double[3] { positionX, positionY, positionZ };
            double coef = CalcCoef(position0, position1);
            double[] directionVector = CalculateDirectionVector(position0, position1, coef);

            double tb = CalcFirstSegmentTime(coef, V, t1);
            Console.WriteLine(tb);
            if (tb > t1 * 0.5) return false;
            double acc = CalculateLinearAcceleration(tb, V);
            if (acc <= 0) return false;

            double[] firstEquation = CalculateParametersOfFirstEquation(acc);
            double[] secondEquation = CalculateParametersOfSecondEquation(V, coef, t1);
            double[] thirdEquation = CalculateParametersOfThirdEquation(coef, acc, t1);


            Console.WriteLine();
            Console.WriteLine("Results: ");
            Console.WriteLine();
            Console.WriteLine("The equation of the first segment of movement is U(t) =" + ShowValue(firstEquation[0]) + "*t^2 + " + ShowValue(firstEquation[1]) + "*t + "+ShowValue(firstEquation[2]));
            Console.WriteLine("The equation of the second segment of movement is U(t) = " + ShowValue(secondEquation[0]) + " * t ^ 2 + " + ShowValue(secondEquation[1]) + " * t + "+ ShowValue(secondEquation[2]));
            Console.WriteLine("The equation of the third segment of movement is U(t) = " + ShowValue(thirdEquation[0]) + " * t ^ 2 + " + ShowValue(thirdEquation[1]) + " * t + "+ ShowValue(thirdEquation[2]));
            Console.WriteLine("The time tb equals "+tb+"s");
            Console.WriteLine();


            if (runSecondPart)
            {
                Console.WriteLine();
                double deltaT = GetInputs("Input the value of delta t, for the discreet representation of the trajector (2.2) (in seconds): ", false);

                double[,] stepMatrix = Part2Point2JointCoords(Part2Point2(t1, tb, acc, deltaT, V, directionVector, position0),new double[3] {theta1, theta2, theta3});

                Console.WriteLine();
                Console.WriteLine("Results (1.2): ");
                Console.WriteLine();

                Console.Write("The list of consecutive positions for joint 1 ([theta1(t0), theta1(t0+deltat)+...+theta1(t1)]): [" + ShowValue(stepMatrix[0,0]));
                for (int i = 1; i < stepMatrix.GetLength(0); i++)
                {
                    Console.Write(", " + ShowValue(stepMatrix[i,0]));
                }
                Console.WriteLine("]");

                Console.Write("The list of consecutive positions for joint 2 ([theta2(t0), theta2(t0+deltat)+...+theta2(t1)]): [" + ShowValue(stepMatrix[0, 1]));
                for (int i = 1; i < stepMatrix.GetLength(0); i++)
                {
                    Console.Write(", " + ShowValue(stepMatrix[i, 1]));
                }
                Console.WriteLine("]");

                Console.Write("The list of consecutive positions for joint 3 ([theta3(t0), theta3(t0+deltat)+...+theta3(t1)]): [" + ShowValue(stepMatrix[0, 2]));
                for (int i = 1; i < stepMatrix.GetLength(0); i++)
                {
                    Console.Write(", " + ShowValue(stepMatrix[i, 2]));
                }
                Console.WriteLine("]");
            }




            return true;
        }




        public static double[,] Part2Point2( double t1, double tb, double acc, double deltat,double Vmax, double[] k, double[] p0)
        {

            int m = 1;
            double[] result = new double[(int)Math.Ceiling(t1 / deltat) + 1];
            result[0] = (acc * deltat * deltat) * 0.5;
            int mSwitch=0;
            while (m * deltat <= t1)
            {
                if (m * deltat <= tb)
                {
                    result[m] = result[m - 1] + (m * deltat * acc) + (acc * deltat * deltat * 0.5);
                }
                if(m* deltat > tb && m * deltat <= t1 - tb)
                {
                    result[m] = result[m - 1] + Vmax * deltat;
                }
                if (m * deltat > t1 -tb)
                {
                    if (mSwitch == 0) mSwitch=m-1;
                    result[m] = result[m - 1] + (Vmax - (m-mSwitch) * deltat * acc) - (acc * deltat * deltat * 0.5);
                }
                m++;
            }

            double[,] carthesian = Part2Point2CarthesianCoords(result, k, p0);

            return carthesian;
        }


        public static double[,] Part2Point2CarthesianCoords(double[] U, double[] k, double[] p0)
        {
            double[,] result = new double[U.Length,3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < U.Length; j++)
                {
                    result[j, i] = p0[i] + U[j] * k[i];

                    Console.WriteLine(result[j,i]);
                }
            }
            return result;
        }

        public static double[,] Part2Point2JointCoords(double[,] carthesian, double[] thetas0)
        {
            double[,] result = new double[carthesian.GetLength(0), carthesian.GetLength(1)];
            result[0, 0] = thetas0[0];
            result[0, 1] = thetas0[1];
            result[0, 2] = thetas0[2];

            for (int i = 1; i < result.GetLength(0); i++)
            {
                List<double[]> stepConfigs = Part1CalcConfigs(result[i - 1, 0], carthesian[i, 0], carthesian[i, 1], carthesian[i, 2]);
                double[] best = SelectCorrectOption(result[i - 1, 0], result[i - 1, 1], result[i - 1, 2], stepConfigs);
                result[i, 0] = best[0];
                result[i, 1] = best[1];
                result[i, 2] = best[2];
            }
            return result;
        }
    }
}
