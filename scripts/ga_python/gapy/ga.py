import time

import numpy as np

try:
    import psutil
except ImportError:  # pragma: no cover - psutil is optional for the GA core.
    psutil = None


def _get_option(gaoptions, key, default=None):
    """
    Lê uma opção de configuração do algoritmo genético.

    Aceita tanto:
    - dicionário: gaoptions["PopulationSize"]
    - objeto com atributos: gaoptions.PopulationSize

    Parameters
    ----------
    gaoptions : dict ou objeto
        Estrutura contendo as opções do GA.
    key : str
        Nome da opção desejada.
    default : any, optional
        Valor a ser retornado caso a opção não exista.

    Returns
    -------
    any
        Valor da opção, se existir; caso contrário, o valor default.
    """
    if gaoptions is None:
        return default

    if isinstance(gaoptions, dict):
        return gaoptions.get(key, default)

    return getattr(gaoptions, key, default)


def randipow(xmax, xpow, n, m=1, rng=None):
    """
    Gera inteiros aleatórios com viés para índices menores.

    Contexto dentro do GA
    ---------------------
    Após ordenar a população por fitness, os melhores indivíduos ficam
    nas primeiras posições (índices menores). Esta função sorteia índices
    com maior probabilidade de escolher esses melhores indivíduos, sem
    eliminar completamente a chance de escolher outros.

    Se xpow > 1:
        a distribuição fica mais concentrada perto de 0.

    Parameters
    ----------
    xmax : int
        Limite superior exclusivo dos índices.
        Os índices gerados estarão em [0, xmax-1].
    xpow : float
        Potência usada para enviesar a distribuição.
    n : int
        Número de linhas da matriz de saída.
    m : int, optional
        Número de colunas da matriz de saída.
    rng : numpy.random.Generator, optional
        Gerador de números aleatórios.

    Returns
    -------
    numpy.ndarray
        Matriz de inteiros com shape (n, m).
    """
    if xmax <= 0:
        raise ValueError("xmax must be > 0.")

    if rng is None:
        rng = np.random.default_rng()

    random_values = rng.random((n, m))
    biased_values = random_values ** xpow
    indices = np.floor(xmax * biased_values).astype(np.int64)

    return indices


def gago(ffit, nbits, gaoptions=None, *args):
    """
    Algoritmo Genético binário didático.

    Convenção adotada
    -----------------
    Este algoritmo assume um problema de MINIMIZAÇÃO:
        fitness menor = indivíduo melhor

    Representação
    -------------
    Cada indivíduo é um vetor binário de tamanho `nbits`.

    Etapas por geração
    ------------------
    1. Avaliação de fitness
    2. Ordenação da população
    3. Seleção de pais
    4. Crossover
    5. Mutação
    6. Elitismo
    7. Repetição por várias gerações

    Parameters
    ----------
    ffit : callable
        Função de fitness.
        Recebe um vetor binário 1D (um indivíduo) e deve retornar
        um valor escalar numérico.
    nbits : int
        Número de bits por indivíduo.
    gaoptions : dict ou objeto, optional
        Configurações do algoritmo genético.

        Opções suportadas
        -----------------
        PopulationSize : int
            Tamanho da população.
        Generations : int
            Número de gerações.
        InitialPopulation : array-like de shape (n_individuos, nbits)
            População inicial. Se não for fornecida, é gerada aleatoriamente.
        MutationRate : float
            Probabilidade de mutação por gene (bit), entre 0 e 1.
        MutationFcn : float
            Alias aceito por compatibilidade. Se MutationRate não existir,
            o código usa MutationFcn.
        EliteCount : int
            Quantidade de melhores indivíduos preservados a cada geração.
        SelectionPressure : float
            Intensidade do viés de seleção. Valores > 1 favorecem mais
            fortemente os melhores indivíduos.
        CrossoverRate : float
            Probabilidade de aplicar crossover entre dois pais.
            Se não ocorrer crossover, o filho é cópia de um dos pais.
        RandomState : int
            Semente aleatória para reprodutibilidade.
        Verbose : bool
            Se True, imprime o melhor fitness por geração.

    *args :
        Mantido por compatibilidade de assinatura, embora não seja usado aqui.

    Returns
    -------
    best_individual : numpy.ndarray
        Melhor indivíduo encontrado (vetor binário 1D).
    final_population : numpy.ndarray
        População final ordenada por fitness crescente.
    final_fitness : numpy.ndarray
        Fitness da população final, em ordem crescente.
    history : dict
        Histórico da execução contendo:
        - "best_fitness_history"
        - "mean_fitness_history"
        - "best_individual_history"
    """
    # ============================================================
    # 1) Leitura e validação básica dos parâmetros
    # ============================================================
    if gaoptions is None:
        gaoptions = {}

    n_bits = int(nbits)
    if n_bits <= 0:
        raise ValueError("nbits must be a positive integer.")

    population_size = _get_option(gaoptions, "PopulationSize", 100)
    n_generations = _get_option(gaoptions, "Generations", 300)

    # Aceita MutationRate como nome principal.
    # Se não existir, tenta MutationFcn por compatibilidade com seu código antigo.
    mutation_rate = _get_option(gaoptions, "MutationRate", None)
    if mutation_rate is None:
        mutation_rate = _get_option(gaoptions, "MutationFcn", 0.01)

    elite_count = _get_option(gaoptions, "EliteCount", 2)
    initial_population = _get_option(gaoptions, "InitialPopulation", None)
    selection_pressure = _get_option(gaoptions, "SelectionPressure", 1.3)
    crossover_rate = _get_option(gaoptions, "CrossoverRate", 1.0)
    random_state = _get_option(gaoptions, "RandomState", None)
    verbose = _get_option(gaoptions, "Verbose", False)

    population_size = int(population_size)
    n_generations = int(n_generations)
    elite_count = int(elite_count)
    selection_pressure = float(selection_pressure)
    crossover_rate = float(crossover_rate)
    mutation_rate = float(mutation_rate)

    if population_size <= 0:
        raise ValueError("PopulationSize must be a positive integer.")

    if n_generations < 0:
        raise ValueError("Generations must be >= 0.")

    if elite_count < 0:
        raise ValueError("EliteCount must be >= 0.")

    if elite_count > population_size:
        elite_count = population_size

    if not (0.0 <= mutation_rate <= 1.0):
        raise ValueError("MutationRate must be in [0, 1].")

    if not (0.0 <= crossover_rate <= 1.0):
        raise ValueError("CrossoverRate must be in [0, 1].")

    if selection_pressure <= 0:
        raise ValueError("SelectionPressure must be > 0.")

    rng = np.random.default_rng(random_state)

    # ============================================================
    # 2) Criação / validação da população inicial
    # ============================================================
    if initial_population is None or (
        hasattr(initial_population, "__len__") and len(initial_population) == 0
    ):
        # Geração aleatória de indivíduos binários
        population = rng.integers(
            low=0,
            high=2,
            size=(population_size, n_bits),
            dtype=np.uint8
        )
    else:
        population = np.asarray(initial_population)

        if population.ndim != 2:
            raise ValueError("InitialPopulation must be a 2D array.")

        if population.shape[1] != n_bits:
            raise ValueError(
                f"InitialPopulation must have exactly {n_bits} columns."
            )

        if population.shape[0] != population_size:
            raise ValueError(
                "PopulationSize does not match InitialPopulation row count."
            )

        # Garante representação binária 0/1
        population = (population != 0).astype(np.uint8)

    # ============================================================
    # 3) Estruturas para armazenar fitness e histórico
    # ============================================================

    best_fitness_history = np.empty(n_generations, dtype=float)
    mean_fitness_history = np.empty(n_generations, dtype=float)
    best_individual_history = np.empty((n_generations, n_bits), dtype=np.uint8)
    generation_elapsed_seconds_history = np.empty(n_generations, dtype=float)
    cumulative_elapsed_seconds_history = np.empty(n_generations, dtype=float)
    process_cpu_seconds_delta_history = np.empty(n_generations, dtype=float)
    process_cpu_seconds_cumulative_history = np.empty(n_generations, dtype=float)
    process_cpu_percent_total_capacity_history = np.empty(n_generations, dtype=float)
    process_cpu_core_percent_history = np.empty(n_generations, dtype=float)
    process_rss_mb_history = np.full(n_generations, np.nan, dtype=float)
    system_memory_percent_history = np.full(n_generations, np.nan, dtype=float)

    logical_cpu_count = psutil.cpu_count(logical=True) if psutil is not None else None
    if not logical_cpu_count:
        logical_cpu_count = 1

    process = psutil.Process() if psutil is not None else None
    ga_started_at = time.perf_counter()
    previous_generation_finished_at = ga_started_at
    process_cpu_started_at = time.process_time()
    previous_process_cpu_seconds = process_cpu_started_at

    # ============================================================
    # 4) Função interna para avaliar a população
    # ============================================================
    def evaluate_population(current_population):
        """
        Avalia todos os indivíduos da população com a função de fitness.
        """
        def evaluate_individual(individual):
            fitness = ffit(individual)

            if not np.isscalar(fitness):
                raise ValueError("ffit must return a scalar fitness value.")

            return fitness

        return np.fromiter(
            (evaluate_individual(individual) for individual in current_population),
            dtype=float,
            count=current_population.shape[0],
        )

    # ============================================================
    # 5) Loop evolutivo principal
    # ============================================================
    for generation_index in range(n_generations):
        # --------------------------------------------------------
        # 5.1 Avalia a população atual
        # --------------------------------------------------------
        fitness_values = evaluate_population(population)

        # --------------------------------------------------------
        # 5.2 Ordena do melhor para o pior
        # Como é minimização, menor fitness vem primeiro.
        # --------------------------------------------------------
        sorted_indices = np.argsort(fitness_values)
        population = population[sorted_indices, :]
        fitness_values = fitness_values[sorted_indices]

        # --------------------------------------------------------
        # 5.3 Salva histórico
        # --------------------------------------------------------
        mean_fitness = float(np.mean(fitness_values))
        best_fitness_history[generation_index] = fitness_values[0]
        mean_fitness_history[generation_index] = mean_fitness
        best_individual_history[generation_index, :] = population[0, :]

        if verbose:
            print(
                f"Geração {generation_index + 1:4d} | "
                f"Melhor fitness = {fitness_values[0]:.10f} | "
                f"Média fitness = {mean_fitness:.10f}"
            )

        # --------------------------------------------------------
        # 5.4 Cria nova população
        # --------------------------------------------------------
        new_population = np.empty_like(population)

        # --------------------------------------------------------
        # 5.5 Elitismo
        # Copia os melhores indivíduos diretamente
        # --------------------------------------------------------
        if elite_count > 0:
            new_population[:elite_count, :] = population[:elite_count, :]

        # --------------------------------------------------------
        # 5.6 Geração dos demais indivíduos por reprodução
        # --------------------------------------------------------
        offspring_count = population_size - elite_count
        if offspring_count > 0:
            # Seleção de dois pais com viés para melhores indivíduos
            parent_indices = randipow(
                xmax=population_size,
                xpow=selection_pressure,
                n=offspring_count,
                m=2,
                rng=rng
            )

            if n_bits > 1:
                crossover_draws = rng.random(offspring_count)
                crossover_points = rng.integers(1, n_bits, size=offspring_count)
                fallback_parent1_draws = rng.random(offspring_count)

                for offspring_index in range(offspring_count):
                    new_index = elite_count + offspring_index
                    parent1 = population[parent_indices[offspring_index, 0], :]
                    parent2 = population[parent_indices[offspring_index, 1], :]

                    if crossover_draws[offspring_index] < crossover_rate:
                        crossover_point = crossover_points[offspring_index]
                        new_population[new_index, :crossover_point] = parent1[:crossover_point]
                        new_population[new_index, crossover_point:] = parent2[crossover_point:]
                    elif fallback_parent1_draws[offspring_index] < 0.5:
                        new_population[new_index, :] = parent1
                    else:
                        new_population[new_index, :] = parent2
            else:
                fallback_parent1_draws = rng.random(offspring_count)

                for offspring_index in range(offspring_count):
                    new_index = elite_count + offspring_index
                    selected_parent_column = 0 if fallback_parent1_draws[offspring_index] < 0.5 else 1
                    new_population[new_index, :] = population[
                        parent_indices[offspring_index, selected_parent_column],
                        :,
                    ]

            # ----------------------------------------------------
            # Mutação por gene
            # Cada bit tem probabilidade mutation_rate de inverter
            # ----------------------------------------------------
            mutation_mask = rng.random((offspring_count, n_bits)) < mutation_rate
            offspring = new_population[elite_count:, :]
            offspring[mutation_mask] = 1 - offspring[mutation_mask]

        # Atualiza população
        population = new_population

        generation_finished_at = time.perf_counter()
        process_cpu_seconds = time.process_time()
        generation_elapsed_seconds = generation_finished_at - previous_generation_finished_at
        process_cpu_seconds_delta = process_cpu_seconds - previous_process_cpu_seconds

        generation_elapsed_seconds_history[generation_index] = generation_elapsed_seconds
        cumulative_elapsed_seconds_history[generation_index] = generation_finished_at - ga_started_at
        process_cpu_seconds_delta_history[generation_index] = process_cpu_seconds_delta
        process_cpu_seconds_cumulative_history[generation_index] = (
            process_cpu_seconds - process_cpu_started_at
        )

        if generation_elapsed_seconds > 0.0:
            process_cpu_core_percent_history[generation_index] = (
                100.0 * process_cpu_seconds_delta / generation_elapsed_seconds
            )
            process_cpu_percent_total_capacity_history[generation_index] = (
                process_cpu_core_percent_history[generation_index] / logical_cpu_count
            )
        else:
            process_cpu_core_percent_history[generation_index] = 0.0
            process_cpu_percent_total_capacity_history[generation_index] = 0.0

        if process is not None:
            process_rss_mb_history[generation_index] = process.memory_info().rss / (1024 ** 2)
            system_memory_percent_history[generation_index] = psutil.virtual_memory().percent

        previous_generation_finished_at = generation_finished_at
        previous_process_cpu_seconds = process_cpu_seconds

    # ============================================================
    # 6) Avaliação final da última população
    # ============================================================
    final_evaluation_started_at = time.perf_counter()
    final_evaluation_process_cpu_started_at = time.process_time()
    fitness_values = evaluate_population(population)
    final_evaluation_finished_at = time.perf_counter()
    final_evaluation_process_cpu_finished_at = time.process_time()
    final_evaluation_elapsed_seconds = (
        final_evaluation_finished_at - final_evaluation_started_at
    )
    final_evaluation_process_cpu_seconds_delta = (
        final_evaluation_process_cpu_finished_at - final_evaluation_process_cpu_started_at
    )

    sorted_indices = np.argsort(fitness_values)
    population = population[sorted_indices, :]
    fitness_values = fitness_values[sorted_indices]

    best_individual = population[0, :].copy()

    history = {
        "best_fitness_history": best_fitness_history,
        "mean_fitness_history": mean_fitness_history,
        "best_individual_history": best_individual_history,
        "generation_elapsed_seconds": generation_elapsed_seconds_history,
        "cumulative_elapsed_seconds": cumulative_elapsed_seconds_history,
        "process_cpu_seconds_delta": process_cpu_seconds_delta_history,
        "process_cpu_seconds_cumulative": process_cpu_seconds_cumulative_history,
        "process_cpu_percent_total_capacity": process_cpu_percent_total_capacity_history,
        "process_cpu_core_percent": process_cpu_core_percent_history,
        "process_rss_mb": process_rss_mb_history,
        "system_memory_percent": system_memory_percent_history,
        "final_evaluation_elapsed_seconds": final_evaluation_elapsed_seconds,
        "final_evaluation_process_cpu_seconds_delta": (
            final_evaluation_process_cpu_seconds_delta
        ),
        "final_process_rss_mb": (
            process.memory_info().rss / (1024 ** 2) if process is not None else np.nan
        ),
        "telemetry_metadata": {
            "logical_cpu_count": int(logical_cpu_count),
            "wall_clock_timer": "time.perf_counter",
            "process_cpu_timer": "time.process_time",
            "ram_source": "psutil.Process.memory_info().rss" if psutil is not None else None,
            "system_memory_source": "psutil.virtual_memory().percent" if psutil is not None else None,
        },
    }

    return best_individual.astype(np.uint8), population.astype(np.uint8), fitness_values, history
