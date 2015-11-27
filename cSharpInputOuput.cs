using System;
using System.IO;
using System.Numerics;
using System.Linq;
namespace Codeforces
{
    class Program
    {
        static void Main(string[] args)
        {
            InputReader cin = new InputReader();
            OutputWriter cout = new OutputWriter();
            //read file and write file
            //InputReader cin = new InputReader(new System.IO.StreamReader(@"C:\Users\yiqiwu\Desktop\Code\log.txt"));
            //OutputWriter cout = new OutputWriter(new System.IO.StreamWriter(@"C:\Users\yiqiwu\Desktop\Code\log.txt"));
            while(cin.HasNext())
            {
                int a = cin.NextInt();
                int b = cin.NextInt();
                cout.WriteLine(a + b);
            }
            cout.Close();
        }
    }
    #region ReaderClass
    class InputReader
    {
        private readonly TextReader reader;
        private string[] buffer = new string[0];
        private int i;
        public InputReader()
        {
            reader = Console.In; 
        }
        public InputReader(TextReader reader)
        {
            this.reader = reader;
            buffer = new string[0];
        }
        public void Close()
        {
            reader.Close();
        }
        public Boolean ReadLine()
        {
            try
            {
                buffer = reader.ReadLine().Split(new[] {' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                i = 0;
            }
            catch (Exception e)
            {
                return false;
            }
            return true;
        }
        public Boolean HasNext()
        {
            while(buffer.Length == 0 || i == buffer.Length)
            {
                if(!ReadLine())
                {
                    return false;
                }
            }
            return true;
        }
        public string NextString()
        {
            while(buffer.Length == 0 || i == buffer.Length)
            {
                ReadLine();
            }
            return buffer[i++];
        }
        public double NextDouble()
        {
            return double.Parse(NextString());
        }

        public int NextInt()
        {
            return int.Parse(NextString());
        }

        public long NextLong()
        {
            return long.Parse(NextString());
        }
        public BigInteger NextBigInteger()
        {
            return BigInteger.Parse(NextString());
        }
        public int[] NextIntArray()
        {
            ReadLine();
            return buffer.Select(int.Parse).ToArray();
        }
        public double[] NextDoubleArray()
        {
            ReadLine();
            return buffer.Select(double.Parse).ToArray();
        }
        public long[] NextLongArray()
        {
            ReadLine();
            return buffer.Select(long.Parse).ToArray();
        }
        public string[] NextStringArray()
        {
            ReadLine();
            return buffer;
        }
    }
    #endregion
    #region WriterClass
    class OutputWriter
    {
        TextWriter writer;
        public OutputWriter()
        {
            writer = Console.Out;
        }
        public OutputWriter(TextWriter writer)
        {
            this.writer = writer;
        }
        public void WriteLine(object a)
        {
            writer.WriteLine(a.ToString());
        }
        public void Write(object a)
        {
            writer.WriteLine(a.ToString());
        }
        public void Close()
        {
            writer.Close();
        }
    }
    #endregion
}
