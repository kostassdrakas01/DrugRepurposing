from typing import Dict, List, Set
import networkx as nx

class BooleanNetworkSimulator:
    def __init__(self, graph: nx.DiGraph):
        self.graph = graph
        self.states = {} # node -> bool

    def initialize_states(self, initial_values: Dict[str, bool] = None, seed: int = 42):
        import random
        random.seed(seed)
        for node in self.graph.nodes():
            if initial_values and node in initial_values:
                self.states[node] = initial_values[node]
            else:
                # Random initial state for simulation
                self.states[node] = random.random() > 0.5

    def update(self, forced_states: Dict[str, bool] = None):
        new_states = {}
        for node in self.graph.nodes():
            if forced_states and node in forced_states:
                new_states[node] = forced_states[node]
                continue

            in_edges = self.graph.in_edges(node, data=True)
            if not in_edges:
                new_states[node] = self.states[node]
                continue

            activators = [src for src, dst, data in in_edges if data.get('weight', 1.0) > 0]
            inhibitors = [src for src, dst, data in in_edges if data.get('weight', 1.0) < 0]

            is_active = any(self.states[act] for act in activators) if activators else self.states[node]
            is_inhibited = any(self.states[inh] for inh in inhibitors)

            new_states[node] = is_active and not is_inhibited

        changed = new_states != self.states
        self.states = new_states
        return changed

    def simulate(self, forced_states: Dict[str, bool] = None, steps: int = 10):
        for _ in range(steps):
            if not self.update(forced_states=forced_states):
                break
        return self.states

    def predict_shift(self, drug_targets: Dict[str, float]):
        """
        Predicts steady state shift.
        drug_targets: node_id -> impact (-1.0 for inhibition, 1.0 for activation)
        """
        # 1. Baseline
        self.initialize_states()
        baseline_state = self.simulate(steps=20).copy()

        # 2. Perturbed State
        forced = {node: (impact > 0) for node, impact in drug_targets.items()}
        # For inhibition, we force it to False. For activation, force it to True.

        self.initialize_states()
        drug_state = self.simulate(forced_states=forced, steps=20).copy()

        shifted_nodes = []
        for node in drug_state:
            if drug_state[node] != baseline_state[node]:
                shifted_nodes.append({
                    "node": node,
                    "from": baseline_state[node],
                    "to": drug_state[node]
                })

        return shifted_nodes
