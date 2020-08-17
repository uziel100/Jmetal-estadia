package org.uma.jmetal.runner.multiobjective;

import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.algorithm.multiobjective.moead.AbstractMOEAD;
import org.uma.jmetal.algorithm.multiobjective.moead.MOEADBuilder;
import org.uma.jmetal.algorithm.multiobjective.moead.MOEADBuilder.Variant;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.impl.crossover.DifferentialEvolutionCrossover;
import org.uma.jmetal.operator.impl.mutation.PolynomialMutation;
import org.uma.jmetal.problem.DoubleProblem;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.AbstractAlgorithmRunner;
import org.uma.jmetal.util.AlgorithmRunner;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.ProblemUtils;

import java.io.FileNotFoundException;
import java.util.List;

/**
 * Class for configuring and running the MOEA/D algorithm
 *
 * @author Antonio J. Nebro <antonio@lcc.uma.es>
 */
public class MOEADConstraintRunner extends AbstractAlgorithmRunner {

    /**
     * @param args Command line arguments.
     * @throws SecurityException Invoking command: java
     * org.uma.jmetal.runner.multiobjective.MOEADRunner problemName
     * [referenceFront]
     */
    public static void main(String[] args) throws FileNotFoundException {
        DoubleProblem problem;
        Algorithm<List<DoubleSolution>> algorithm;
        MutationOperator<DoubleSolution> mutation;
        DifferentialEvolutionCrossover crossover;

        String problemName;
        String referenceParetoFront = "";
        if (args.length == 1) {
            problemName = args[0];
        } else if (args.length == 2) {
            problemName = args[0];
            referenceParetoFront = args[1];
        } else {
            // Binh2
            // ConstrEx
            // Osyczka2
            // Srinivas
            // Tanaka
            // Golinski
            problemName = "org.uma.jmetal.problem.multiobjective.Golinski";
            //referenceParetoFront = "jmetal-problem/src/test/resources/pareto_fronts/Tanaka.pf";
            referenceParetoFront = "C:\\Users\\UZIEL\\Documents\\Estadia-JMETAL\\JMetal-master\\src\\pareto_fronts\\Golinski.pf";
        }

        problem = (DoubleProblem) ProblemUtils.<DoubleSolution>loadProblem(problemName);

        double cr = 0.9;
        double f = 0.9;
        crossover = new DifferentialEvolutionCrossover(cr, f, "rand/1/bin");

        double mutationProbability = 1.0 / problem.getNumberOfVariables();
        double mutationDistributionIndex = 20.0;
        mutation = new PolynomialMutation(mutationProbability, mutationDistributionIndex);

        int tam = 30;
        String[][] outputMatriz = new String[tam][3];
        for (int i = 0; i < tam; i++) {

            algorithm = new MOEADBuilder(problem, Variant.ConstraintMOEAD)
                    .setCrossover(crossover)
                    .setMutation(mutation)
                    .setMaxEvaluations(1500)
                    .setPopulationSize(100)
                    .setResultPopulationSize(100)
                    .setNeighborhoodSelectionProbability(0.9)
                    .setMaximumNumberOfReplacedSolutions(2)
                    .setNeighborSize(20)
                    .setFunctionType(AbstractMOEAD.FunctionType.TCHE)
                    .setDataDirectory("MOEAD_Weights")
                    .build();

            AlgorithmRunner algorithmRunner = new AlgorithmRunner.Executor(algorithm)
                    .execute();

            List<DoubleSolution> population = algorithm.getResult();
            long computingTime = algorithmRunner.getComputingTime();

            //JMetalLogger.logger.info("Total execution time: " + computingTime + "ms");
            printFinalSolutionSet(population);
            if (!referenceParetoFront.equals("")) {
                //printQualityIndicators(population, referenceParetoFront);

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

        System.out.println("---------------------------------------------------");
        System.out.println("IDG");

        prom = 0;

        for (int i = 0; i < tam; i++) {
            prom += Double.parseDouble(outputMatriz[i][2]);
            //System.out.println(outputMatriz[i][2]);
        }

        System.out.println(prom / tam);

    }
}
