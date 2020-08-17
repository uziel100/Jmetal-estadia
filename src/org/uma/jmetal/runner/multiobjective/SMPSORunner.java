package org.uma.jmetal.runner.multiobjective;

import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.algorithm.multiobjective.smpso.SMPSOBuilder;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.impl.mutation.PolynomialMutation;
import org.uma.jmetal.problem.DoubleProblem;
import org.uma.jmetal.problem.multiobjective.*;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.AbstractAlgorithmRunner;
import org.uma.jmetal.util.AlgorithmRunner;
import org.uma.jmetal.util.archive.BoundedArchive;
import org.uma.jmetal.util.archive.impl.CrowdingDistanceArchive;
import org.uma.jmetal.util.evaluator.impl.SequentialSolutionListEvaluator;

import java.util.List;

/**
 * Class for configuring and running the SMPSO algorithm
 *
 * @author Antonio J. Nebro <antonio@lcc.uma.es>
 */
public class SMPSORunner extends AbstractAlgorithmRunner {

    /**
     * @param args Command line arguments. The first (optional) argument
     * specifies the problem to solve.
     * @throws org.uma.jmetal.util.JMetalException
     * @throws java.io.IOException
     * @throws SecurityException Invoking command: java
     * org.uma.jmetal.runner.multiobjective.SMPSORunner problemName
     * [referenceFront]
     */
    public static void main(String[] args) throws Exception {
        DoubleProblem problem;
        Algorithm<List<DoubleSolution>> algorithm;
        MutationOperator<DoubleSolution> mutation;

        String referenceParetoFront = "";

        String problemName;
        if (args.length == 1) {
            problemName = args[0];
        } else if (args.length == 2) {
            problemName = args[0];
            referenceParetoFront = args[1];
        } else {

            // Cambiar la Ruta de base
            //String base = "C:\\Users\\UZIEL\\Documents\\Estadia-JMETAL";            
            //referenceParetoFront = base + "\\JMetal-master\\src\\pareto_fronts\\ZDT1.pf";
            problemName = "org.uma.jmetal.problem.multiobjective.Tanaka";
            //referenceParetoFront = "jmetal-problem/src/test/resources/pareto_fronts/Tanaka.pf";
            referenceParetoFront = "C:\\Users\\UZIEL\\Documents\\Estadia-JMETAL\\JMetal-master\\src\\pareto_fronts\\Tanaka.pf";
        }
        //System.out.println("Warning: the problem name is not used anymore and may be removed later.") ;
        System.out.println("Warning: current problem name: " + problemName);

        problem = new Tanaka();

        BoundedArchive<DoubleSolution> archive = new CrowdingDistanceArchive<DoubleSolution>(100);

        double mutationProbability = 1.0 / problem.getNumberOfVariables();
        double mutationDistributionIndex = 20.0;
        mutation = new PolynomialMutation(mutationProbability, mutationDistributionIndex);

        int tam = 30;
        String[][] outputMatriz = new String[tam][3];

        for (int i = 0; i < tam; i++) {
            algorithm = new SMPSOBuilder(problem, archive)
                    .setMutation(mutation)
                    .setMaxIterations(15)
                    .setSwarmSize(100)
                    .setSolutionListEvaluator(new SequentialSolutionListEvaluator<DoubleSolution>())
                    .build();

            AlgorithmRunner algorithmRunner = new AlgorithmRunner.Executor(algorithm)
                    .execute();

            List<DoubleSolution> population = algorithm.getResult();
            long computingTime = algorithmRunner.getComputingTime();

            //JMetalLogger.logger.info("Total execution time: " + computingTime + "ms");
            printFinalSolutionSet(population);
            if (!referenceParetoFront.equals("")) {
                String[] metricas = printQualityIndicators(population, referenceParetoFront);
                outputMatriz[i][0] = metricas[0];
                outputMatriz[i][1] = metricas[1];
                outputMatriz[i][2] = metricas[2];
            }
        }

        System.out.println("Hipevolumen");
        double prom = 0;
        for (int i = 0; i < tam; i++) {
            prom += Double.parseDouble(outputMatriz[i][0]);
            //System.out.println(outputMatriz[i][0]);
        }

        System.out.println(prom / tam);

        System.out.println("---------------------------------------------------");
        System.out.println("Epsilon");
        prom = 0;
        for (int i = 0; i < tam; i++) {
            prom += Double.parseDouble(outputMatriz[i][1]);
            //System.out.println(outputMatriz[i][1]);
        }

        System.out.println(prom / tam);

        System.out.println("---------------------------------------------------*");
        System.out.println("IDG");

        prom = 0;

        for (int i = 0; i < tam; i++) {
            prom += Double.parseDouble(outputMatriz[i][2]);
            //System.out.println(outputMatriz[i][2]);
        }

        System.out.println(prom / tam);

    }
}
