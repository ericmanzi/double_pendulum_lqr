//
// Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
//
// File:      lo1_fusion.cs
//
// Purpose: Demonstrates how to solve the problem
//
// minimize 3*x1 + 5*x2 + x3 + x3
// such that
//          3*x1 + 2*x2        +   x3 = 30,
//          2*x1 + 3*x2 +   x3 +   x3 > 15,
//                      + 3*x3 + 2*x3 < 25
// and
//          x0,x1,x2 > 0,
//          0 < x3 < 10
//
using System;
using mosek.fusion;
namespace mosek
{
  namespace fusion
  {
    namespace example
    {
      public class lo1
      {
        public static void Main(string[] args)
        {
          double[][] A =
            { new double[] { 3.0, 2.0, 0.0, 1.0 },
              new double[] { 2.0, 3.0, 1.0, 1.0 },
              new double[] { 0.0, 0.0, 3.0, 2.0 } };
          double[] c = { 3.0, 5.0, 1.0, 1.0  };

          // Create a model with the name 'lo1'
          Model M = new Model("lo1");
            
          // Create variable 'x' of length 4
          Variable x = M.Variable("x", 4, Domain.GreaterThan(0.0));
            
          // Create three constraints
          M.Constraint("c1", Expr.Dot(A[0], x), Domain.EqualsTo(30.0));
          M.Constraint("c2", Expr.Dot(A[1], x), Domain.GreaterThan(15.0));
          M.Constraint("c3", Expr.Dot(A[2], x), Domain.LessThan(25.0));
            
          // Set the objective function to (c^t * x)
          M.Objective("obj", ObjectiveSense.Maximize, Expr.Dot(c, x));
                  
          // Solve the problem
          M.Solve();
            
          // Get the solution values
          double[] sol = x.Level();
          Console.WriteLine("[x0,x1,x2,x3] = [{0}, {1}, {2}, {3} ]",sol[0],sol[1],sol[2],sol[3]);
        }
      }
    }
  }
}
   
