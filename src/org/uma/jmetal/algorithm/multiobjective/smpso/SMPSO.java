package org.uma.jmetal.algorithm.multiobjective.smpso;

import org.uma.jmetal.algorithm.impl.AbstractParticleSwarmOptimization;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.problem.DoubleProblem;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.archive.BoundedArchive;
import org.uma.jmetal.util.comparator.DominanceComparator;
import org.uma.jmetal.util.evaluator.SolutionListEvaluator;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;
import org.uma.jmetal.util.solutionattribute.impl.GenericSolutionAttribute;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import org.uma.jmetal.solution.impl.ArrayDoubleSolution;
import sun.security.krb5.internal.crypto.Des3;

/**
 * This class implements the SMPSO algorithm described in: SMPSO: A new
 * PSO-based metaheuristic for multi-objective optimization MCDM 2009. DOI:
 * http://dx.doi.org/10.1109/MCDM.2009.4938830
 *
 * @author Antonio J. Nebro <antonio@lcc.uma.es>
 */
@SuppressWarnings("serial")
public class SMPSO extends AbstractParticleSwarmOptimization<DoubleSolution, List<DoubleSolution>> {

    private DoubleProblem problem;

    private double c1Max;
    private double c1Min;
    private double c2Max;
    private double c2Min;
    private double r1Max;
    private double r1Min;
    private double r2Max;
    private double r2Min;
    private double weightMax;
    private double weightMin;
    private double changeVelocity1;
    private double changeVelocity2;

    private int swarmSize;
    private int maxIterations;
    private int iterations;

    private GenericSolutionAttribute<DoubleSolution, DoubleSolution> localBest;
    private double[][] speed;

    private JMetalRandom randomGenerator;

    private BoundedArchive<DoubleSolution> leaders;
    private Comparator<DoubleSolution> dominanceComparator;

    private MutationOperator<DoubleSolution> mutation;

    private double deltaMax[];
    private double deltaMin[];

    private SolutionListEvaluator<DoubleSolution> evaluator;

    private enum _MethodForPosition {
        Evolutionay, Reflection, Random, Wrapping, Boundary, Conservatism, CentroidV1, Infinity
    }

    private enum _MethodForVelocity {
        None, AbsorbZero, DeterministBack, RandomBack, Adjust
    }

    private final _MethodForPosition _position = _MethodForPosition.CentroidV1;
    private final _MethodForVelocity _velocity = _MethodForVelocity.AbsorbZero;

    /**
     * Constructor
     */
    public SMPSO(DoubleProblem problem, int swarmSize, BoundedArchive<DoubleSolution> leaders,
            MutationOperator<DoubleSolution> mutationOperator, int maxIterations, double r1Min, double r1Max,
            double r2Min, double r2Max, double c1Min, double c1Max, double c2Min, double c2Max,
            double weightMin, double weightMax, double changeVelocity1, double changeVelocity2,
            SolutionListEvaluator<DoubleSolution> evaluator) {
        this.problem = problem;
        this.swarmSize = swarmSize;
        this.leaders = leaders;
        this.mutation = mutationOperator;
        this.maxIterations = maxIterations;

        this.r1Max = r1Max;
        this.r1Min = r1Min;
        this.r2Max = r2Max;
        this.r2Min = r2Min;
        this.c1Max = c1Max;
        this.c1Min = c1Min;
        this.c2Max = c2Max;
        this.c2Min = c2Min;
        this.weightMax = weightMax;
        this.weightMin = weightMin;
        this.changeVelocity1 = changeVelocity1;
        this.changeVelocity2 = changeVelocity2;

        randomGenerator = JMetalRandom.getInstance();
        this.evaluator = evaluator;

        dominanceComparator = new DominanceComparator<DoubleSolution>();
        localBest = new GenericSolutionAttribute<DoubleSolution, DoubleSolution>();
        speed = new double[swarmSize][problem.getNumberOfVariables()];

        deltaMax = new double[problem.getNumberOfVariables()];
        deltaMin = new double[problem.getNumberOfVariables()];
        for (int i = 0; i < problem.getNumberOfVariables(); i++) {
            deltaMax[i] = (problem.getUpperBound(i) - problem.getLowerBound(i)) / 2.0;
            deltaMin[i] = -deltaMax[i];
        }
    }

    protected void updateLeadersDensityEstimator() {
        leaders.computeDensityEstimator();
    }

    @Override
    protected void initProgress() {
        iterations = 1;
        updateLeadersDensityEstimator();
    }

    @Override
    protected void updateProgress() {
        iterations += 1;
        updateLeadersDensityEstimator();
    }

    @Override
    protected boolean isStoppingConditionReached() {
        return iterations >= maxIterations;
    }

    @Override
    protected List<DoubleSolution> createInitialSwarm() {
        List<DoubleSolution> swarm = new ArrayList<>(swarmSize);

        DoubleSolution newSolution;
        for (int i = 0; i < swarmSize; i++) {
            newSolution = problem.createSolution();
            swarm.add(newSolution);
        }

        return swarm;
    }

    @Override
    protected List<DoubleSolution> evaluateSwarm(List<DoubleSolution> swarm) {
        swarm = evaluator.evaluate(swarm, problem);

        return swarm;
    }

    @Override
    protected void initializeLeader(List<DoubleSolution> swarm) {
        for (DoubleSolution particle : swarm) {
            leaders.add(particle);
        }
    }

    @Override
    protected void initializeVelocity(List<DoubleSolution> swarm) {
        for (int i = 0; i < swarm.size(); i++) {
            for (int j = 0; j < problem.getNumberOfVariables(); j++) {
                speed[i][j] = 0.0;
            }
        }
    }

    @Override
    protected void initializeParticlesMemory(List<DoubleSolution> swarm) {
        for (DoubleSolution particle : swarm) {
            localBest.setAttribute(particle, (DoubleSolution) particle.copy());
        }
    }

    @Override
    protected void updateVelocity(List<DoubleSolution> swarm) {
        double r1, r2, c1, c2;
        double wmax, wmin;
        DoubleSolution bestGlobal;

        for (int i = 0; i < swarm.size(); i++) {
            DoubleSolution particle = (DoubleSolution) swarm.get(i).copy();
            DoubleSolution bestParticle = (DoubleSolution) localBest.getAttribute(swarm.get(i)).copy();

            bestGlobal = selectGlobalBest();

            r1 = randomGenerator.nextDouble(r1Min, r1Max);
            r2 = randomGenerator.nextDouble(r2Min, r2Max);
            c1 = randomGenerator.nextDouble(c1Min, c1Max);
            c2 = randomGenerator.nextDouble(c2Min, c2Max);
            wmax = weightMax;
            wmin = weightMin;

            for (int var = 0; var < particle.getNumberOfVariables(); var++) {
                speed[i][var] = velocityConstriction(constrictionCoefficient(c1, c2) * (inertiaWeight(iterations, maxIterations, wmax, wmin) * speed[i][var]
                        + c1 * r1 * (bestParticle.getVariableValue(var) - particle.getVariableValue(var))
                        + c2 * r2 * (bestGlobal.getVariableValue(var) - particle.getVariableValue(var))),
                        deltaMax, deltaMin, var);
            }
        }
    }
    
    
    @Override
    protected void updatePosition(List<DoubleSolution> swarm) {        
        double _particleRepaired = 0.0;
        double _velocityRepaired = 0.0;    
    
        for (int i = 0; i < swarmSize; i++) {            
            DoubleSolution particle = swarm.get(i);
            DoubleSolution particleBefore = (DoubleSolution) particle.copy();

            for (int j = 0; j < particle.getNumberOfVariables(); j++) {
                particle.setVariableValue(j, particle.getVariableValue(j) + speed[i][j]);
                
                
                if (particle.getVariableValue(j) < problem.getLowerBound(j)) {
                    particle = _handlePosition(particleBefore, particle, j, problem.getLowerBound(j));
                    _velocityRepaired = _handleVelocity(particleBefore, particle, j, speed[i][j]);
                    /*_particleRepaired = _handlePosition(particleBefore, particle, j, problem.getLowerBound(j));
                    _velocityRepaired = _handleVelocity(particleBefore, particle, j, speed[i][j]);
                    */
                    speed[i][j] = _velocityRepaired;
                    //particle.setVariableValue(j, _particleRepaired);                    
                }
                if (particle.getVariableValue(j) > problem.getUpperBound(j)) {
                    particle = _handlePosition(particleBefore, particle, j, problem.getUpperBound(j));
                    _velocityRepaired = _handleVelocity(particleBefore, particle, j, speed[i][j]);
                     
                    /*_particleRepaired = _handlePosition(particleBefore, particle, j, problem.getUpperBound(j));
                    _velocityRepaired = _handleVelocity(particleBefore, particle, j, speed[i][j]);*/

                    speed[i][j] = _velocityRepaired;
                    //particle.setVariableValue(j, _particleRepaired);                    
                }
                
            }
        }        
    } 

    //METODOS PARA POSICION
    
    private DoubleSolution _handlePosition(DoubleSolution particleBefore, DoubleSolution particle, int i, double violatedLimitValue) {
        double particleRepaired = 0.0;
        DoubleSolution particleRep = (DoubleSolution) particle.copy();
        
        switch (this._position) {
            case Boundary:
                particleRepaired = _methodBoundary(violatedLimitValue);
                break;
            case Reflection:
                particleRepaired = _methodReflection(particle.getVariableValue(i), violatedLimitValue);
                break;
            case Random:
                particleRepaired = _methodRandom(particle, i);
                break;
            case Evolutionay:
                particleRepaired = _methodEvolutionay(violatedLimitValue, i);
                break;
            case Wrapping:
                particleRepaired = _methodWrapping(particle, i);
                break;
            case Conservatism:
                particleRepaired = particleBefore.getVariableValue(i);
                break;
            case CentroidV1:
                particleRep = _methodCentroidV1(particle, i);
                break;
        }
        if(this._position != _MethodForPosition.CentroidV1)
            particleRep.setVariableValue(i, particleRepaired);
        
        return particleRep;
    }

    private double _methodBoundary(double violatedLimitValue) {
        return violatedLimitValue;
    }

    private double _methodReflection(double currentValue, double violatedLimitValue) {
        return 2 * violatedLimitValue - currentValue;
    }

    private double _methodRandom(DoubleSolution particle, int i) {
        return particle.getLowerBound(i) + randomGenerator.nextDouble() * (particle.getUpperBound(i) - particle.getLowerBound(i));
    }

    private double _methodEvolutionay(double violatedLimitValue, int i) {
        double rnd = randomGenerator.nextDouble();
        DoubleSolution gbest = selectGlobalBest();

        return rnd * violatedLimitValue + (1 - rnd) * gbest.getVariableValue(i);
    }

    private double _methodWrapping(DoubleSolution particle, int i) {
        double currentValue = particle.getVariableValue(i);
        double lowetLimit = particle.getLowerBound(i);
        double upperLimit = particle.getUpperBound(i);

        if (currentValue < lowetLimit)
            return upperLimit - (lowetLimit - currentValue) % (upperLimit - lowetLimit);        

        return lowetLimit + (currentValue - upperLimit) % (upperLimit - lowetLimit);
    }

    
    private DoubleSolution _methodCentroidV1(DoubleSolution particle, int i){
        DoubleSolution gbest = selectGlobalBest();
        //DoubleSolution particleRepairedDim = (DoubleSolution) particle.copy();
        DoubleSolution particleNew = (DoubleSolution) particle.copy();
        
        double operation = 0;
        double particleRepairedByRandom = _methodRandom(particle, i);
        particle.setVariableValue(i, particleRepairedByRandom);
        
        for (int d = 0; d < particle.getNumberOfVariables(); d++) {
            operation = (particle.getVariableValue(d) + gbest.getVariableValue(d)) / 2;
            particleNew.setVariableValue(d, operation);
        }
        
        return particleNew;    
    }
    
    
    //METODOS PARA LA VELOCIDAD
    private double _handleVelocity(DoubleSolution particleBefore, DoubleSolution particle, int i, double currentVelocity) {
        double velocityRepaired = 0.0;

        switch (this._velocity) {
            case None:
                velocityRepaired = currentVelocity;
                break;
            case AbsorbZero:
                velocityRepaired = 0;
                break;
            case DeterministBack:
                velocityRepaired = _methodDeterministBack(particle, i, currentVelocity);
                break;
            case RandomBack:
                velocityRepaired = _methodRandomBack(currentVelocity);
                break;
            case Adjust:
                velocityRepaired = _methodAdjust(particleBefore, particle, i);
                break;
        }

        return velocityRepaired;
    }

    private double _methodDeterministBack(DoubleSolution particle, int i, double currentVelocity) {
        double value = (particle.getVariableValue(i) < particle.getLowerBound(i)) ? this.changeVelocity1 : this.changeVelocity2;

        return value * currentVelocity;
    }

    private double _methodRandomBack(double currentVelocity) {
        return -randomGenerator.nextDouble() * currentVelocity;
    }

    private double _methodAdjust(DoubleSolution particleBefore, DoubleSolution currentParticle, int i) {
        return  currentParticle.getVariableValue(i) - particleBefore.getVariableValue(i);
    }

    @Override
    protected void perturbation(List<DoubleSolution> swarm) {
        for (int i = 0; i < swarm.size(); i++) {
            if ((i % 6) == 0) {
                mutation.execute(swarm.get(i));
            }
        }
    }

    @Override
    protected void updateLeaders(List<DoubleSolution> swarm) {
        for (DoubleSolution particle : swarm) {
            leaders.add((DoubleSolution) particle.copy());
        }
    }

    @Override
    protected void updateParticlesMemory(List<DoubleSolution> swarm) {
        for (int i = 0; i < swarm.size(); i++) {
            int flag = dominanceComparator.compare(swarm.get(i), localBest.getAttribute(swarm.get(i)));
            if (flag != 1) {
                DoubleSolution particle = (DoubleSolution) swarm.get(i).copy();
                localBest.setAttribute(swarm.get(i), particle);
            }
        }
    }

    @Override
    public List<DoubleSolution> getResult() {
        return leaders.getSolutionList();
    }

    protected DoubleSolution selectGlobalBest() {
        DoubleSolution one, two;
        DoubleSolution bestGlobal;
        int pos1 = randomGenerator.nextInt(0, leaders.getSolutionList().size() - 1);
        int pos2 = randomGenerator.nextInt(0, leaders.getSolutionList().size() - 1);
        one = leaders.getSolutionList().get(pos1);
        two = leaders.getSolutionList().get(pos2);

        if (leaders.getComparator().compare(one, two) < 1) {
            bestGlobal = (DoubleSolution) one.copy();
        } else {
            bestGlobal = (DoubleSolution) two.copy();
        }

        return bestGlobal;
    }

    private double velocityConstriction(double v, double[] deltaMax, double[] deltaMin,
            int variableIndex) {
        double result;

        double dmax = deltaMax[variableIndex];
        double dmin = deltaMin[variableIndex];

        result = v;

        if (v > dmax) {
            result = dmax;
        }

        if (v < dmin) {
            result = dmin;
        }

        return result;
    }

    protected double constrictionCoefficient(double c1, double c2) {
        double rho = c1 + c2;
        if (rho <= 4) {
            return 1.0;
        } else {
            return 2 / (2 - rho - Math.sqrt(Math.pow(rho, 2.0) - 4.0 * rho));
        }
    }

    private double inertiaWeight(int iter, int miter, double wma, double wmin) {
        return wma;
    }

    @Override
    public String getName() {
        return "SMPSO";
    }

    @Override
    public String getDescription() {
        return "Speed contrained Multiobjective PSO";
    }

    /* Getters */
    public int getSwarmSize() {
        return swarmSize;
    }

    public int getMaxIterations() {
        return maxIterations;
    }

    public int getIterations() {
        return iterations;
    }

    /* Setters */
    public void setIterations(int iterations) {
        this.iterations = iterations;
    }
}
