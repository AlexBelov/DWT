using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Drawing;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;

namespace Wavelet
{
    class Program
    {
        class DWT
        {
            private Bitmap inImage;
            private double[,] matrix;
            private double[,] outMatrix;
            private double[,] beautifiedMatrix;
            //private double[] CL = { 1 / Math.Sqrt(2), 1 / Math.Sqrt(2) };
            private double[] CL = { 
                (1 + Math.Sqrt(3)) / (4 * Math.Sqrt(2)),
                (3 + Math.Sqrt(3)) / (4 * Math.Sqrt(2)),
                (3 - Math.Sqrt(3)) / (4 * Math.Sqrt(2)),
                (1 - Math.Sqrt(3)) / (4 * Math.Sqrt(2)) 
            };
            private double[] CH;

            public DWT(string path)
            {
                inImage = new Bitmap(path);
                matrix = new double[inImage.Width, inImage.Height];
                outMatrix = new double[inImage.Width, inImage.Height];
                beautifiedMatrix = new double[inImage.Width, inImage.Height];
                BitmapToMatrix();
                HpfCoeffs();
            }

            private void HpfCoeffs()
            {
                int N = CL.Length;
                CH = new double[N];
                for (int i = 0; i < N; i++)
                {
                    CH[i] = Math.Pow(-1, i) * CL[N - i - 1];
                }
            }

            private void BitmapToMatrix()
            {
                for (int x = 0; x < inImage.Width; x++)
                {
                    for (int y = 0; y < inImage.Height; y++)
                    {
                        Color C = inImage.GetPixel(x, y);
                        int Value = (C.R + C.G + C.B) / 3;
                        matrix[x, y] = Value / 255.0;
                    }
                }
            }

            private double[] Pconv(double[] data, int delta = 0)
            {
                double[] result = new double[data.Length];
                int N = CL.Length;
                int M = data.Length;
                int iResult = 0;

                for (int k = 0; k < M; k += 2)
                {
                    double sL = 0;
                    double sH = 0;

                    for (int i = 0; i < N; i++)
                    {
                        sL += data[(k + i - delta) % M] * CL[i];
                        sH += data[(k + i - delta) % M] * CH[i];
                    }

                    result[iResult++] = sL;
                    result[iResult++] = sH;
                }

                return result;
            }

            private void CalculateDWT()
            {
                double[,] matrixT = new double[inImage.Width, inImage.Height];
                Array.Copy(matrix, matrixT, inImage.Width * inImage.Height);

                // Rows
                for (int i = 0; i < inImage.Height; i++)
                {
                    double[] row = new double[inImage.Width];
                    for (int k = 0; k < inImage.Width; k++)
                    {
                        row[k] = matrixT[i, k];
                    }

                    double[] pconvRow = new double[inImage.Width];
                    pconvRow = Pconv(row);

                    for (int k = 0; k < inImage.Width; k++)
                    {
                        matrixT[i, k] = pconvRow[k];
                    }
                }

                // Columns
                for (int i = 0; i < inImage.Width; i++)
                {
                    double[] column = new double[inImage.Height];
                    for (int k = 0; k < inImage.Height; k++)
                    {
                        column[k] = matrixT[k, i];
                    }

                    double[] pconvColumn = new double[inImage.Height];
                    pconvColumn = Pconv(column);

                    for (int k = 0; k < inImage.Height; k++)
                    {
                        matrixT[k, i] = pconvColumn[k];
                    }
                }

                Array.Copy(matrixT, outMatrix, inImage.Width * inImage.Height);
            }

            private void Beautify()
            {
                int dim = inImage.Width;
                double[,] data = new double[dim, dim];

                int N = dim * dim / 4;
                int temp_i = 0;
                double[] temp_pixels = new double[N];

                temp_i = 0;
                temp_pixels = new double[N];
                for (int i = 0; i < dim; i += 2)
                {
                    for (int j = 0; j < dim; j += 2)
                    {
                        temp_pixels[temp_i++] = outMatrix[i, j];
                    }
                }
                temp_i = 0;
                for (int i = 0; i < dim / 2; i += 1)
                {
                    for (int j = 0; j < dim / 2; j += 1)
                    {
                        data[i, j] = temp_pixels[temp_i++];
                    }
                }

                temp_i = 0;
                temp_pixels = new double[N];
                for (int i = 1; i < dim; i += 2)
                {
                    for (int j = 0; j < dim; j += 2)
                    {
                        temp_pixels[temp_i++] = outMatrix[i, j];
                    }
                }
                temp_i = 0;
                for (int i = dim / 2; i < dim; i += 1)
                {
                    for (int j = 0; j < dim / 2; j += 1)
                    {
                        data[i, j] = temp_pixels[temp_i++];
                    }
                }

                temp_i = 0;
                temp_pixels = new double[N];
                for (int i = 0; i < dim; i += 2)
                {
                    for (int j = 1; j < dim; j += 2)
                    {
                        temp_pixels[temp_i++] = outMatrix[i, j];
                    }
                }
                temp_i = 0;
                for (int i = 0; i < dim / 2; i += 1)
                {
                    for (int j = dim / 2; j < dim; j += 1)
                    {
                        data[i, j] = temp_pixels[temp_i++];
                    }
                }

                temp_i = 0;
                temp_pixels = new double[N];
                for (int i = 1; i < dim; i += 2)
                {
                    for (int j = 1; j < dim; j += 2)
                    {
                        temp_pixels[temp_i++] = outMatrix[i, j];
                    }
                }
                temp_i = 0;
                for (int i = dim / 2; i < dim; i += 1)
                {
                    for (int j = dim / 2; j < dim; j += 1)
                    {
                        data[i, j] = temp_pixels[temp_i++];
                    }
                }

                Array.Copy(data, beautifiedMatrix, inImage.Width * inImage.Height);
            }

            public void Processing()
            {
                CalculateDWT();
                Beautify();
            }

            public void SaveImage()
            {
                Bitmap outImage = new Bitmap(inImage.Width, inImage.Height);

                for (int x = 0; x < outImage.Width; x++)
                {
                    for (int y = 0; y < outImage.Height; y++)
                    {
                        double outputPixel = Math.Abs(beautifiedMatrix[x, y] * 140);
                        if (outputPixel > 255)
                        {
                            outputPixel = 255;
                        }
                        int e = System.Convert.ToInt32(outputPixel);
                        Color pixel = Color.FromArgb(e, e, e);
                        outImage.SetPixel(x, y, pixel);
                    }
                }

                outImage.Save("Intermediate.png");
            }

            public void SaveMatrix()
            {
                Serialize(outMatrix, "Intermediate.dwt");
            }

            private static void Serialize(object t, string path)
            {
                using (Stream stream = File.Open(path, FileMode.Create))
                {
                    BinaryFormatter bformatter = new BinaryFormatter();
                    bformatter.Serialize(stream, t);
                }
            }
        }

        class IDWT
        {
            private int dim;
            private double[,] matrix;
            private double[,] outMatrix;
            //private double[] CL = { 1 / Math.Sqrt(2), 1 / Math.Sqrt(2) };
            private double[] CL = { 
                (1 + Math.Sqrt(3)) / (4 * Math.Sqrt(2)),
                (3 + Math.Sqrt(3)) / (4 * Math.Sqrt(2)),
                (3 - Math.Sqrt(3)) / (4 * Math.Sqrt(2)),
                (1 - Math.Sqrt(3)) / (4 * Math.Sqrt(2)) 
            };
            private double[] CH;
            private double[] iCL;
            private double[] iCH;

            public IDWT(string path)
            {
                matrix = (double[,])Deserialize(path);
                HpfCoeffs();
                ICoeffs();
                dim = matrix.GetLength(0);
                outMatrix = new double[dim, dim];
            }

            private void HpfCoeffs()
            {
                int N = CL.Length;
                CH = new double[N];
                for (int i = 0; i < N; i++)
                {
                    CH[i] = Math.Pow(-1, i) * CL[N - i - 1];
                }
            }

            private void ICoeffs()
            {
                int N = CL.Length;
                iCH = new double[N];
                iCL = new double[N];
                int i = 0;
                for (int k = 0; k < CL.Length; k += 2)
                {
                    int index_1 = k - 2;
                    if (index_1 < 0)
                    {
                        index_1 = CL.Length - k - 2;
                    }
                    int index_2 = k - 1;
                    if (index_2 < 0)
                    {
                        index_2 = CL.Length - k - 1;
                    }
                    iCL[i] = CL[index_1];
                    iCH[i] = CL[index_2];
                    i++;
                    iCL[i] = CH[index_1];
                    iCH[i] = CH[index_2];
                    i++;
                }
            }

            private static object Deserialize(string path)
            {
                using (Stream stream = File.Open(path, FileMode.Open))
                {
                    BinaryFormatter bformatter = new BinaryFormatter();
                    return bformatter.Deserialize(stream);
                }
            }

            private double[] Pconv(double[] data, int delta = 0)
            {
                double[] result = new double[data.Length];
                int N = iCL.Length;
                int M = data.Length;
                int iResult = 0;

                for (int k = 0; k < M; k += 2)
                {
                    double sL = 0;
                    double sH = 0;

                    for (int i = 0; i < N; i++)
                    {
                        int index = k + i - delta;
                        if (index < 0)
                        {
                            index = data.Length - k - i + delta;
                        }
                        sL += data[index % M] * iCL[i];
                        sH += data[index % M] * iCH[i];
                    }

                    result[iResult++] = sL;
                    result[iResult++] = sH;
                }

                return result;
            }

            public void Processing()
            {
                CalculateIDWT();
            }

            private void CalculateIDWT()
            {
                double[,] matrixT = new double[dim, dim];
                Array.Copy(matrix, matrixT, dim * dim);

                // Columns
                for (int i = 0; i < dim; i++)
                {
                    double[] column = new double[dim];
                    for (int k = 0; k < dim; k++)
                    {
                        column[k] = matrixT[k, i];
                    }

                    double[] pconvColumn = new double[dim];
                    pconvColumn = Pconv(column, iCL.Length - 2);

                    for (int k = 0; k < dim; k++)
                    {
                        matrixT[k, i] = pconvColumn[k];
                    }
                }

                // Rows
                for (int i = 0; i < dim; i++)
                {
                    double[] row = new double[dim];
                    for (int k = 0; k < dim; k++)
                    {
                        row[k] = matrixT[i, k];
                    }

                    double[] pconvRow = new double[dim];
                    pconvRow = Pconv(row, iCL.Length - 2);

                    for (int k = 0; k < dim; k++)
                    {
                        matrixT[i, k] = pconvRow[k];
                    }
                }

                Array.Copy(matrixT, outMatrix, dim * dim);
            }

            public void SaveImage()
            {
                Bitmap outImage = new Bitmap(dim, dim);

                for (int x = 0; x < outImage.Width; x++)
                {
                    for (int y = 0; y < outImage.Height; y++)
                    {
                        double outputPixel = Math.Abs(outMatrix[x, y] * 255);
                        if (outputPixel > 255)
                        {
                            outputPixel = 255;
                        }
                        int e = System.Convert.ToInt32(outputPixel);
                        Color pixel = Color.FromArgb(e, e, e);
                        outImage.SetPixel(x, y, pixel);
                    }
                }

                outImage.Save("Output.png");
            }
        }

        static void Main(string[] args)
        {
            //DWT dwt = new DWT("Boat.png");
            DWT dwt = new DWT("Lenna.jpg");
            dwt.Processing();
            dwt.SaveImage();
            dwt.SaveMatrix();

            IDWT iDwt = new IDWT("Intermediate.dwt");
            iDwt.Processing();
            iDwt.SaveImage();

            //Console.WriteLine("\n" + "Press any key to close...");
            //Console.ReadKey();
        }
    }
}
