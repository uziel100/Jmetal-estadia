package org.uma.jmetal.operator.impl.crossover;

import org.uma.jmetal.operator.CrossoverOperator;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.JMetalLogger;
import org.uma.jmetal.util.pseudorandom.BoundedRandomGenerator;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;
import org.uma.jmetal.util.pseudorandom.RandomGenerator;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import org.uma.jmetal.algorithm.multiobjective.moead.util.MOEADUtils;
import org.uma.jmetal.util.comparator.impl.ViolationThresholdComparator;
import org.uma.jmetal.util.solutionattribute.impl.NumberOfViolatedConstraints;
import org.uma.jmetal.util.solutionattribute.impl.OverallConstraintViolation;

/**
 * Differential evolution crossover operator
 *
 * @author Antonio J. Nebro
 *
 * Comments: - The operator receives two parameters: the current individual and
 * an array of three parent individuals - The best and rand variants depends on
 * the third parent, according whether it represents the current of the "best"
 * individual or a random one. The implementation of both variants are the same,
 * due to that the parent selection is external to the crossover operator. -
 * Implemented variants: - rand/1/bin (best/1/bin) - rand/1/exp (best/1/exp) -
 * current-to-rand/1 (current-to-best/1) - current-to-rand/1/bin
 * (current-to-best/1/bin) - current-to-rand/1/exp (current-to-best/1/exp)
 */
@SuppressWarnings("serial")
public class DifferentialEvolutionCrossover implements CrossoverOperator<DoubleSolution> {

    private static final double DEFAULT_CR = 0.5;
    private static final double DEFAULT_F = 0.5;
    private static final double DEFAULT_K = 0.5;
    private static final String DEFAULT_DE_VARIANT = "rand/1/bin";

    private double cr;
    private double f;
    private double k;
    // DE variant (rand/1/bin, rand/1/exp, etc.)
    private String variant;

    private DoubleSolution currentSolution;

    private BoundedRandomGenerator<Integer> jRandomGenerator;
    private BoundedRandomGenerator<Double> crRandomGenerator;

    private JMetalRandom randomGenerator;
    private OverallConstraintViolation<DoubleSolution> overallConstraintViolation;
    private NumberOfViolatedConstraints<DoubleSolution> numberOfViolatedConstraints;
    private ViolationThresholdComparator<DoubleSolution> violationThresholdComparator;
    private List<DoubleSolution> populationCopy;

    /**
     * Constructor
     */
    public DifferentialEvolutionCrossover() {
        this(DEFAULT_CR, DEFAULT_F, DEFAULT_K, DEFAULT_DE_VARIANT);
    }

    /**
     * Constructor
     *
     * @param cr
     * @param f
     * @param variant
     */
    public DifferentialEvolutionCrossover(double cr, double f, String variant) {
        this(cr, f, variant, (a, b) -> JMetalRandom.getInstance().nextInt(a, b), (a, b) -> JMetalRandom.getInstance().nextDouble(a, b));
        randomGenerator = JMetalRandom.getInstance();

        overallConstraintViolation = new OverallConstraintViolation<DoubleSolution>();
        numberOfViolatedConstraints = new NumberOfViolatedConstraints<DoubleSolution>();
        violationThresholdComparator = new ViolationThresholdComparator<>();
    }

    /**
     * Constructor
     *
     * @param cr
     * @param f
     * @param variant
     * @param jRandomGenerator
     * @param crRandomGenerator
     */
    public DifferentialEvolutionCrossover(double cr, double f, String variant, RandomGenerator<Double> randomGenerator) {
        this(cr, f, variant, BoundedRandomGenerator.fromDoubleToInteger(randomGenerator), BoundedRandomGenerator.bound(randomGenerator));
    }

    /**
     * Constructor
     *
     * @param cr
     * @param f
     * @param variant
     * @param jRandomGenerator
     * @param crRandomGenerator
     */
    public DifferentialEvolutionCrossover(double cr, double f, String variant, BoundedRandomGenerator<Integer> jRandomGenerator, BoundedRandomGenerator<Double> crRandomGenerator) {
        this.cr = cr;        
        this.f = f;
        this.k = DEFAULT_K;
        this.variant = variant;

        this.jRandomGenerator = jRandomGenerator;
        this.crRandomGenerator = crRandomGenerator;
    }

    /**
     * Constructor
     */
    public DifferentialEvolutionCrossover(double cr, double f, double k, String variant) {
        this(cr, f, variant);
        this.k = k;
    }

    /* Getters */
    public double getCr() {
        return cr;
    }

    public double getF() {
        return f;
    }

    public double getK() {
        return k;
    }

    public String getVariant() {
        return variant;
    }

    /* Setters */
    public void setCurrentSolution(DoubleSolution current) {
        this.currentSolution = current;
    }

    public void setCr(double cr) {
        this.cr = cr;
    }

    public void setF(double f) {
        this.f = f;
    }

    public void setK(double k) {
        this.k = k;
    }

    //*************************************
    private ArrayList<DoubleSolution> SFS; // conjunto de soluciones factibles
    private ArrayList<DoubleSolution> SIS; // conjunto de soluciones no factibles

    private int AFS; // cantidad de soluciones factibles

    private int[] probability;    
    private boolean particleRapaired;
    private int methodUsed;
    

    public void _setSFSandSIS(List<DoubleSolution> population) {
        SFS = new ArrayList<>();
        SIS = new ArrayList<>();
        AFS = 0;

        for (int i = 0; i < population.size(); i++) {

            if (!isFeasible(population.get(i))) {
                SIS.add((DoubleSolution) population.get(i).copy());
            } else {
                SFS.add((DoubleSolution) population.get(i).copy());
                AFS++;
            }
        }
    }

    private boolean isFeasible(DoubleSolution solution) {
        return (numberOfViolatedConstraints.getAttribute(solution) == 0);
    }

    public void _setPopulation(List<DoubleSolution> population) {
        this.populationCopy = population;
    }
    
    public void _setProbability(int[] arr){
        this.probability = arr;
    }
    
    public int _getMethodUsed(){
        return methodUsed;
    }
    
    public boolean _isParticleRepaired(){
        return this.particleRapaired;
    }


    //*************************************
    /**
     * Execute() method
     * @param parentSolutions
     * @return 
     */
    @Override
    public List<DoubleSolution> execute(List<DoubleSolution> parentSolutions) {
        
        DoubleSolution child = (DoubleSolution) currentSolution.copy();
        particleRapaired = false;         

        // STEP 4. Checking the DE variant
        if ((DEFAULT_DE_VARIANT.equals(variant))
                || "best/1/bin".equals(variant)) {

            //genera el vector mutante y aplica el operador de cruza            
            _setOperationCruza(child, parentSolutions);

            if (!isValid(child)) {                  
                //_methodBoundary(child);
                //_methodReflection(child);
                //_methodEvolutionay(child, parentSolutions);
                //_methodWrapping(child);
                //_methodRandom(child);
                //_methodCentroidV1(child, 2);
                //_methodResAndRand(child);                                
                
                
                _selectMethodForRapaired();                               
                _applyMethod(child, parentSolutions);                                
                particleRapaired = true;
            }

        } else {
            JMetalLogger.logger.severe("DifferentialEvolutionCrossover.execute: "
                    + " unknown DE variant (" + variant + ")");
            Class<String> cls = String.class;
            String name = cls.getName();
            throw new JMetalException("Exception in " + name + ".execute()");
        }

        List<DoubleSolution> result = new ArrayList<>(1);
        result.add(child);
        return result;
    }

    private void _selectMethodForRapaired() {
        if (AFS == 0) {
            methodUsed = 1;
        } else {
            methodUsed = probability[randomGenerator.nextInt(0, 99)];
        }
    }
       

    private void _applyMethod(DoubleSolution child, List<DoubleSolution> parentSolutions) {
        switch (methodUsed) {
            case 1:
                _methodResAndRand(child);                
                //_methodRandom(child);
                break;
            case 2:
                _methodCentroidV1(child, 2);
                break;
            case 3:
                //_methodReflection(child);
                //_methodEvolutionay(child, parentSolutions);
                _methodWrapping(child);
                break;
            case 4:
                //_methodBoundary(child);
                _methodEvolutionay(child, parentSolutions);
                break;
            default:
                break;
        }
    }

  

    private void _methodResAndRand(DoubleSolution child) {
        int nres = 3 * child.getNumberOfVariables();
        int i = 0;

        do {
            _setOperationCruza(child, _getSolutions());
            i++;
        } while ( !isValid(child) || nres > i);

        if (!isValid(child)) {
            _methodRandom(child);
        }

    }

    private List<DoubleSolution> _getSolutions() {
        List<DoubleSolution> solutions = new LinkedList<>();

        int[] permutation = new int[populationCopy.size()];
        MOEADUtils.randomPermutation(permutation, populationCopy.size());

        int rnd = randomGenerator.nextInt(0, populationCopy.size() - 3);

        solutions.add(populationCopy.get(permutation[rnd]));
        solutions.add(populationCopy.get(permutation[rnd + 1]));
        solutions.add(populationCopy.get(permutation[rnd + 2]));

        return solutions;
    }

    private void _setOperationCruza(DoubleSolution child, List<DoubleSolution> parentSolutions) {
        int jrand = jRandomGenerator.getRandomValue(0, child.getNumberOfVariables() - 1);        
        for (int j = 0; j < child.getNumberOfVariables(); j++) {
            if (crRandomGenerator.getRandomValue(0.0, 1.0) < cr || j == jrand) {
                double value;
                value = parentSolutions.get(2).getVariableValue(j) + f * (parentSolutions.get(0).getVariableValue(j)
                        - parentSolutions.get(1).getVariableValue(j));

                child.setVariableValue(j, value);
            } else {
                //cr = 0.0;
                double value;
                value = currentSolution.getVariableValue(j);
                child.setVariableValue(j, value);
            }
        }
    }

    private boolean isValid(DoubleSolution child) {
        for (int j = 0; j < child.getNumberOfVariables(); j++) {
            if (child.getVariableValue(j) < child.getLowerBound(j) || child.getVariableValue(j) > child.getUpperBound(j)) {
                return false;
            }
        }
        return true;
    }

    private void _methodBoundary(DoubleSolution child) {
        for (int j = 0; j < child.getNumberOfVariables(); j++) {
            if (child.getVariableValue(j) < child.getLowerBound(j)) {
                child.setVariableValue(j, child.getLowerBound(j));
            }

            if (child.getVariableValue(j) > child.getUpperBound(j)) {
                child.setVariableValue(j, child.getUpperBound(j));
            }

        }
    }

    private void _methodReflection(DoubleSolution child) {
        double value = 0;
        for (int j = 0; j < child.getNumberOfVariables(); j++) {

            if (child.getVariableValue(j) < child.getLowerBound(j)) {
                value = 2 * child.getLowerBound(j) - child.getVariableValue(j);
                child.setVariableValue(j, value);
            }

            if (child.getVariableValue(j) > child.getUpperBound(j)) {
                value = 2 * child.getUpperBound(j) - child.getVariableValue(j);
                child.setVariableValue(j, value);
            }
        }
    }

    private void _methodRandom(DoubleSolution child) {
        double value = 0;
        for (int j = 0; j < child.getNumberOfVariables(); j++) {
            if (child.getVariableValue(j) < child.getLowerBound(j) || child.getVariableValue(j) > child.getUpperBound(j)) {
                value = child.getLowerBound(j) + randomGenerator.nextDouble() * (child.getUpperBound(j) - child.getLowerBound(j));
                child.setVariableValue(j, value);
            }
        }
    }

    private DoubleSolution _methodRandom2(DoubleSolution child) {
        double value = 0;
        DoubleSolution child2 = (DoubleSolution) child.copy();

        for (int j = 0; j < child2.getNumberOfVariables(); j++) {
            if (child2.getVariableValue(j) < child.getLowerBound(j) || child2.getVariableValue(j) > child2.getUpperBound(j)) {
                value = child.getLowerBound(j) + randomGenerator.nextDouble() * (child2.getUpperBound(j) - child.getLowerBound(j));
                child2.setVariableValue(j, value);
            }
        }
        return child2;
    }

    private void _methodEvolutionay(DoubleSolution child, List<DoubleSolution> parentSolutions) {
        double rnd = 0;
        double value = 0;
        for (int j = 0; j < child.getNumberOfVariables(); j++) {
            rnd = randomGenerator.nextDouble();

            if (child.getVariableValue(j) < child.getLowerBound(j)) {
                value = rnd * child.getLowerBound(j) + (1 - rnd) * parentSolutions.get(2).getVariableValue(j);
                child.setVariableValue(j, value);
            }

            if (child.getVariableValue(j) > child.getUpperBound(j)) {
                value = rnd * child.getUpperBound(j) + (1 - rnd) * parentSolutions.get(2).getVariableValue(j);
                child.setVariableValue(j, value);
            }
        }
    }

    private void _methodWrapping(DoubleSolution child) {
        double value = 0;
        for (int j = 0; j < child.getNumberOfVariables(); j++) {
            if (child.getVariableValue(j) < child.getLowerBound(j)) {
                value = child.getUpperBound(j) - (child.getLowerBound(j) - child.getVariableValue(j)) % (child.getUpperBound(j) - child.getLowerBound(j));
                child.setVariableValue(j, value);
            }

            if (child.getVariableValue(j) > child.getUpperBound(j)) {
                value = child.getLowerBound(j) + (child.getVariableValue(j) - child.getUpperBound(j)) % (child.getUpperBound(j) - child.getLowerBound(j));
                child.setVariableValue(j, value);
            }

        }
    }

    private void _methodCentroidV1(DoubleSolution child, int k) {
        double value = 0;
        DoubleSolution[] childs = new DoubleSolution[k];

        DoubleSolution wp = _getCurrentBestSolution();

        for (int i = 0; i < k; i++) {
            childs[i] = _methodRandom2(child);
        }

        for (int j = 0; j < child.getNumberOfVariables(); j++) {
            //se calcula el centroid
            value = wp.getVariableValue(j);
            for (int i = 0; i < k; i++) {
                value += childs[i].getVariableValue(j);
            }
            value /= (k + 1);
            child.setVariableValue(j, value);
        }
    }

    private DoubleSolution _getCurrentBestSolution() {

        if (AFS > 0 && randomGenerator.nextDouble() > 0.5) {
            return (DoubleSolution) SFS.get(randomGenerator.nextInt(0, SFS.size() - 1)).copy();
        } else {
            if (SIS.isEmpty()) {
                return (DoubleSolution) SFS.get(randomGenerator.nextInt(0, SFS.size() - 1)).copy();
            }

            return _getChildBestSIS();
        }
    }

    public DoubleSolution _getChildBestSIS() {

        DoubleSolution best = (DoubleSolution) SIS.get(0).copy();

        for (int i = 1; i < SIS.size(); i++) {

            if (violationThresholdComparator.compare(SIS.get(i), best) == -1) {
                best = (DoubleSolution) SIS.get(i).copy();
            }
        }

        return best;
    }

    private void _methodCentroidV1(DoubleSolution child, List<DoubleSolution> parentSolutions, int k) {
        double value = 0;
        DoubleSolution[] childs = new DoubleSolution[k];

        for (int i = 0; i < k; i++) {
            childs[i] = _methodRandom2(child);
        }

        for (int j = 0; j < child.getNumberOfVariables(); j++) {
            //se calcula el centroid
            value = parentSolutions.get(2).getVariableValue(j);
            for (int i = 0; i < k; i++) {
                value += childs[i].getVariableValue(j);
            }
            value /= (k + 1);
            child.setVariableValue(j, value);
        }
    }

    @Override
    public int getNumberOfRequiredParents() {
        return 3;
    }

    @Override
    public int getNumberOfGeneratedChildren() {
        return 1;
    }
}
