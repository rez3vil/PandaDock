# -*- coding: utf-8 -*-
"""
Machine Learning-based rescoring module
Implements various ML approaches for pose rescoring and ranking
"""

import numpy as np
from typing import List, Dict, Any, Optional, Tuple
import logging
from pathlib import Path
import pickle
import json

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docking.base_engine import Pose


class MLRescorer:
    """
    Machine Learning-based rescoring system
    
    Features:
    - Multiple ML model support (Random Forest, Neural Network, etc.)
    - Feature extraction from poses
    - Model training and inference
    - Consensus scoring
    - Uncertainty quantification
    """
    
    def __init__(self, config=None):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Model parameters
        self.model_type = 'random_forest'  # Default model
        self.model_path = None
        self.feature_scaler = None
        self.trained_model = None
        
        # Feature extraction parameters
        self.feature_types = [
            'energy_terms',
            'geometric_descriptors',
            'interaction_counts',
            'physicochemical_properties',
            'pocket_descriptors'
        ]
        
        # Model performance metrics
        self.training_metrics = {}
        self.validation_metrics = {}
        
        self.logger.info("Initialized MLRescorer")
    
    def load_model(self, model_path: str):
        """Load pre-trained ML model"""
        self.logger.info(f"Loading ML model from {model_path}")
        
        model_path = Path(model_path)
        
        if not model_path.exists():
            raise FileNotFoundError(f"Model file not found: {model_path}")
        
        try:
            with open(model_path, 'rb') as f:
                model_data = pickle.load(f)
            
            self.trained_model = model_data['model']
            self.feature_scaler = model_data.get('scaler', None)
            self.model_type = model_data.get('model_type', 'unknown')
            self.training_metrics = model_data.get('training_metrics', {})
            
            self.logger.info(f"Successfully loaded {self.model_type} model")
            
        except Exception as e:
            self.logger.error(f"Error loading model: {e}")
            raise
    
    def save_model(self, model_path: str):
        """Save trained ML model"""
        self.logger.info(f"Saving ML model to {model_path}")
        
        if self.trained_model is None:
            raise ValueError("No trained model to save")
        
        model_data = {
            'model': self.trained_model,
            'scaler': self.feature_scaler,
            'model_type': self.model_type,
            'training_metrics': self.training_metrics,
            'validation_metrics': self.validation_metrics
        }
        
        with open(model_path, 'wb') as f:
            pickle.dump(model_data, f)
        
        self.logger.info("Model saved successfully")
    
    def train_model(self, training_data: List[Dict[str, Any]], 
                   validation_data: Optional[List[Dict[str, Any]]] = None):
        """
        Train ML model for pose rescoring
        
        Args:
            training_data: List of training examples with poses and labels
            validation_data: Optional validation data
        """
        self.logger.info(f"Training {self.model_type} model with {len(training_data)} examples")
        
        # Extract features and labels
        X_train, y_train = self._prepare_training_data(training_data)
        
        # Feature scaling
        self.feature_scaler = self._create_feature_scaler()
        X_train_scaled = self.feature_scaler.fit_transform(X_train)
        
        # Train model
        self.trained_model = self._create_model()
        self.trained_model.fit(X_train_scaled, y_train)
        
        # Calculate training metrics
        train_predictions = self.trained_model.predict(X_train_scaled)
        self.training_metrics = self._calculate_metrics(y_train, train_predictions)
        
        # Validation
        if validation_data:
            X_val, y_val = self._prepare_training_data(validation_data)
            X_val_scaled = self.feature_scaler.transform(X_val)
            val_predictions = self.trained_model.predict(X_val_scaled)
            self.validation_metrics = self._calculate_metrics(y_val, val_predictions)
            
            self.logger.info(f"Validation RMSE: {self.validation_metrics['rmse']:.3f}")
            self.logger.info(f"Validation R²: {self.validation_metrics['r2']:.3f}")
        
        self.logger.info(f"Training completed. Training RMSE: {self.training_metrics['rmse']:.3f}")
    
    def _prepare_training_data(self, data: List[Dict[str, Any]]) -> Tuple[np.ndarray, np.ndarray]:
        """Prepare training data by extracting features and labels"""
        features = []
        labels = []
        
        for example in data:
            pose = example['pose']
            label = example['label']  # True binding affinity or score
            
            # Extract features from pose
            pose_features = self.extract_features(pose)
            features.append(pose_features)
            labels.append(label)
        
        return np.array(features), np.array(labels)
    
    def _create_model(self):
        """Create ML model based on model type"""
        if self.model_type == 'random_forest':
            from sklearn.ensemble import RandomForestRegressor
            return RandomForestRegressor(
                n_estimators=100,
                max_depth=10,
                random_state=42,
                n_jobs=-1
            )
        
        elif self.model_type == 'xgboost':
            try:
                import xgboost as xgb
                return xgb.XGBRegressor(
                    n_estimators=100,
                    max_depth=6,
                    learning_rate=0.1,
                    random_state=42
                )
            except ImportError:
                self.logger.warning("XGBoost not available, falling back to Random Forest")
                from sklearn.ensemble import RandomForestRegressor
                return RandomForestRegressor(n_estimators=100, random_state=42)
        
        elif self.model_type == 'neural_network':
            return self._create_neural_network()
        
        else:
            raise ValueError(f"Unknown model type: {self.model_type}")
    
    def _create_neural_network(self):
        """Create neural network model"""
        from sklearn.neural_network import MLPRegressor
        
        return MLPRegressor(
            hidden_layer_sizes=(100, 50),
            activation='relu',
            solver='adam',
            alpha=0.001,
            learning_rate='adaptive',
            max_iter=1000,
            random_state=42
        )
    
    def _create_feature_scaler(self):
        """Create feature scaler"""
        from sklearn.preprocessing import StandardScaler
        return StandardScaler()
    
    def _calculate_metrics(self, y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
        """Calculate performance metrics"""
        from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
        
        rmse = np.sqrt(mean_squared_error(y_true, y_pred))
        mae = mean_absolute_error(y_true, y_pred)
        r2 = r2_score(y_true, y_pred)
        
        # Pearson correlation
        correlation = np.corrcoef(y_true, y_pred)[0, 1]
        
        return {
            'rmse': rmse,
            'mae': mae,
            'r2': r2,
            'correlation': correlation
        }
    
    def rescore_poses(self, poses: List[Pose]) -> List[Pose]:
        """
        Rescore poses using trained ML model
        
        Args:
            poses: List of poses to rescore
            
        Returns:
            List of poses with updated ML scores
        """
        if self.trained_model is None:
            raise ValueError("No trained model available. Load or train a model first.")
        
        self.logger.info(f"Rescoring {len(poses)} poses with ML model")
        
        # Extract features for all poses
        features = []
        for pose in poses:
            pose_features = self.extract_features(pose)
            features.append(pose_features)
        
        # Scale features
        X = np.array(features)
        if self.feature_scaler:
            X_scaled = self.feature_scaler.transform(X)
        else:
            X_scaled = X
        
        # Predict scores
        ml_scores = self.trained_model.predict(X_scaled)
        
        # Update poses with ML scores
        for i, pose in enumerate(poses):
            pose.ml_score = ml_scores[i]
            
            # Optionally replace original score
            if hasattr(self.config, 'use_ml_score_as_primary') and self.config.use_ml_score_as_primary:
                pose.score = ml_scores[i]
        
        # Sort by ML score
        poses.sort(key=lambda x: getattr(x, 'ml_score', x.score))
        
        return poses
    
    def extract_features(self, pose: Pose) -> np.ndarray:
        """
        Extract features from a pose for ML rescoring
        
        Args:
            pose: Pose object
            
        Returns:
            Feature vector
        """
        features = []
        
        # Energy terms
        if 'energy_terms' in self.feature_types:
            energy_features = self._extract_energy_features(pose)
            features.extend(energy_features)
        
        # Geometric descriptors
        if 'geometric_descriptors' in self.feature_types:
            geometric_features = self._extract_geometric_features(pose)
            features.extend(geometric_features)
        
        # Interaction counts
        if 'interaction_counts' in self.feature_types:
            interaction_features = self._extract_interaction_features(pose)
            features.extend(interaction_features)
        
        # Physicochemical properties
        if 'physicochemical_properties' in self.feature_types:
            physicochemical_features = self._extract_physicochemical_features(pose)
            features.extend(physicochemical_features)
        
        # Pocket descriptors
        if 'pocket_descriptors' in self.feature_types:
            pocket_features = self._extract_pocket_features(pose)
            features.extend(pocket_features)
        
        return np.array(features)
    
    def _extract_energy_features(self, pose: Pose) -> List[float]:
        """Extract energy-related features"""
        features = [
            pose.energy,
            pose.vdw_energy,
            pose.electrostatic_energy,
            pose.hbond_energy,
            pose.hydrophobic_energy,
            pose.solvation_energy,
            pose.entropy_energy,
            pose.clash_score
        ]
        
        # Normalize large energy values
        features = [min(max(f, -100.0), 100.0) for f in features]
        
        return features
    
    def _extract_geometric_features(self, pose: Pose) -> List[float]:
        """Extract geometric descriptors"""
        coordinates = pose.coordinates
        
        if len(coordinates) < 2:
            return [0.0] * 10  # Return default features
        
        # Center of mass
        center_of_mass = np.mean(coordinates, axis=0)
        
        # Radius of gyration
        distances_to_center = np.linalg.norm(coordinates - center_of_mass, axis=1)
        radius_of_gyration = np.sqrt(np.mean(distances_to_center**2))
        
        # Molecular dimensions
        min_coords = np.min(coordinates, axis=0)
        max_coords = np.max(coordinates, axis=0)
        dimensions = max_coords - min_coords
        
        # Asphericity (measure of non-spherical shape)
        asphericity = self._calculate_asphericity(coordinates)
        
        # Distance from binding site center (placeholder)
        binding_site_center = np.array([0.0, 0.0, 0.0])  # Would be actual binding site center
        distance_to_binding_site = np.linalg.norm(center_of_mass - binding_site_center)
        
        features = [
            radius_of_gyration,
            dimensions[0],  # x dimension
            dimensions[1],  # y dimension
            dimensions[2],  # z dimension
            asphericity,
            distance_to_binding_site,
            np.std(distances_to_center),  # Spread of atoms
            len(coordinates),  # Number of atoms
            np.mean(distances_to_center),  # Average distance to center
            np.max(distances_to_center)   # Max distance to center
        ]
        
        return features
    
    def _calculate_asphericity(self, coordinates: np.ndarray) -> float:
        """Calculate asphericity (measure of non-spherical shape)"""
        if len(coordinates) < 3:
            return 0.0
        
        # Calculate inertia tensor
        center = np.mean(coordinates, axis=0)
        centered_coords = coordinates - center
        
        inertia_tensor = np.zeros((3, 3))
        for coord in centered_coords:
            r_squared = np.dot(coord, coord)
            inertia_tensor += r_squared * np.eye(3) - np.outer(coord, coord)
        
        # Calculate eigenvalues
        eigenvalues = np.linalg.eigvals(inertia_tensor)
        eigenvalues = np.sort(eigenvalues)[::-1]  # Sort descending
        
        # Asphericity = (λ1 - 0.5*(λ2 + λ3)) / (λ1 + λ2 + λ3)
        if np.sum(eigenvalues) > 0:
            asphericity = (eigenvalues[0] - 0.5 * (eigenvalues[1] + eigenvalues[2])) / np.sum(eigenvalues)
        else:
            asphericity = 0.0
        
        return asphericity
    
    def _extract_interaction_features(self, pose: Pose) -> List[float]:
        """Extract interaction-related features"""
        features = [
            len(pose.hbond_interactions),
            len(pose.hydrophobic_interactions),
            len(pose.salt_bridge_interactions),
            pose.binding_site_coverage,
            pose.ligand_efficiency
        ]
        
        # Average interaction energies
        if pose.hbond_interactions:
            avg_hbond_energy = np.mean([interaction.get('energy', 0.0) for interaction in pose.hbond_interactions])
        else:
            avg_hbond_energy = 0.0
        
        if pose.hydrophobic_interactions:
            avg_hydrophobic_energy = np.mean([interaction.get('energy', 0.0) for interaction in pose.hydrophobic_interactions])
        else:
            avg_hydrophobic_energy = 0.0
        
        features.extend([avg_hbond_energy, avg_hydrophobic_energy])
        
        return features
    
    def _extract_physicochemical_features(self, pose: Pose) -> List[float]:
        """Extract physicochemical property features"""
        # These would be calculated from actual molecular structure
        # For now, using placeholders
        
        features = [
            300.0,  # Molecular weight (placeholder)
            2.5,    # LogP (placeholder)
            3,      # HBD count (placeholder)
            5,      # HBA count (placeholder)
            60.0,   # TPSA (placeholder)
            4       # Rotatable bonds (placeholder)
        ]
        
        return features
    
    def _extract_pocket_features(self, pose: Pose) -> List[float]:
        """Extract binding pocket descriptors"""
        # These would be calculated from actual pocket analysis
        # For now, using placeholders
        
        features = [
            1000.0,  # Pocket volume (placeholder)
            500.0,   # Pocket surface area (placeholder)
            0.7,     # Pocket hydrophobicity (placeholder)
            2.0,     # Pocket depth (placeholder)
            0.3      # Pocket druggability score (placeholder)
        ]
        
        return features
    
    def predict_binding_affinity(self, pose: Pose) -> float:
        """
        Predict binding affinity for a pose
        
        Args:
            pose: Pose object
            
        Returns:
            Predicted binding affinity (kcal/mol)
        """
        if self.trained_model is None:
            raise ValueError("No trained model available")
        
        features = self.extract_features(pose)
        
        if self.feature_scaler:
            features_scaled = self.feature_scaler.transform(features.reshape(1, -1))
        else:
            features_scaled = features.reshape(1, -1)
        
        predicted_affinity = self.trained_model.predict(features_scaled)[0]
        
        return predicted_affinity
    
    def predict_with_uncertainty(self, pose: Pose) -> Tuple[float, float]:
        """
        Predict binding affinity with uncertainty estimation
        
        Args:
            pose: Pose object
            
        Returns:
            Tuple of (predicted_affinity, uncertainty)
        """
        if self.model_type == 'random_forest':
            # Use tree predictions for uncertainty
            features = self.extract_features(pose)
            
            if self.feature_scaler:
                features_scaled = self.feature_scaler.transform(features.reshape(1, -1))
            else:
                features_scaled = features.reshape(1, -1)
            
            # Get predictions from all trees
            tree_predictions = []
            for tree in self.trained_model.estimators_:
                pred = tree.predict(features_scaled)[0]
                tree_predictions.append(pred)
            
            mean_prediction = np.mean(tree_predictions)
            uncertainty = np.std(tree_predictions)
            
            return mean_prediction, uncertainty
        
        else:
            # For other models, use simple prediction
            prediction = self.predict_binding_affinity(pose)
            return prediction, 0.0  # No uncertainty estimation
    
    def consensus_scoring(self, poses: List[Pose], methods: List[str]) -> List[Pose]:
        """
        Apply consensus scoring using multiple methods
        
        Args:
            poses: List of poses
            methods: List of scoring methods to combine
            
        Returns:
            Poses with consensus scores
        """
        self.logger.info(f"Applying consensus scoring with {len(methods)} methods")
        
        for pose in poses:
            scores = []
            
            # Collect scores from different methods
            if 'ml_score' in methods and hasattr(pose, 'ml_score'):
                scores.append(pose.ml_score)
            
            if 'physics_score' in methods:
                scores.append(pose.energy)
            
            if 'vina_score' in methods:
                # Would calculate Vina score
                scores.append(pose.score)
            
            # Calculate consensus score
            if scores:
                pose.consensus_score = np.mean(scores)
            else:
                pose.consensus_score = pose.score
        
        # Sort by consensus score
        poses.sort(key=lambda x: x.consensus_score)
        
        return poses
    
    def get_feature_importance(self) -> Dict[str, float]:
        """Get feature importance from trained model"""
        if self.trained_model is None:
            raise ValueError("No trained model available")
        
        if hasattr(self.trained_model, 'feature_importances_'):
            # For tree-based models
            importances = self.trained_model.feature_importances_
            
            # Create feature names
            feature_names = []
            if 'energy_terms' in self.feature_types:
                feature_names.extend(['energy', 'vdw', 'electrostatic', 'hbond', 'hydrophobic', 'solvation', 'entropy', 'clash'])
            if 'geometric_descriptors' in self.feature_types:
                feature_names.extend(['radius_gyration', 'dim_x', 'dim_y', 'dim_z', 'asphericity', 'dist_binding_site', 'spread', 'num_atoms', 'avg_dist', 'max_dist'])
            if 'interaction_counts' in self.feature_types:
                feature_names.extend(['hbond_count', 'hydrophobic_count', 'salt_bridge_count', 'binding_coverage', 'ligand_efficiency', 'avg_hbond_energy', 'avg_hydrophobic_energy'])
            if 'physicochemical_properties' in self.feature_types:
                feature_names.extend(['molecular_weight', 'logp', 'hbd', 'hba', 'tpsa', 'rotatable_bonds'])
            if 'pocket_descriptors' in self.feature_types:
                feature_names.extend(['pocket_volume', 'pocket_surface', 'pocket_hydrophobicity', 'pocket_depth', 'druggability'])
            
            # Create importance dictionary
            importance_dict = {}
            for i, name in enumerate(feature_names):
                if i < len(importances):
                    importance_dict[name] = importances[i]
            
            return importance_dict
        
        else:
            return {}
    
    def get_rescorer_info(self) -> Dict[str, Any]:
        """Get information about the ML rescorer"""
        info = {
            'model_type': self.model_type,
            'feature_types': self.feature_types,
            'trained': self.trained_model is not None,
            'training_metrics': self.training_metrics,
            'validation_metrics': self.validation_metrics
        }
        
        if self.trained_model is not None:
            info['feature_importance'] = self.get_feature_importance()
        
        return info