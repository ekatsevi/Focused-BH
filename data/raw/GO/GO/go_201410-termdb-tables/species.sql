-- MySQL dump 10.13  Distrib 5.1.61, for unknown-linux-gnu (x86_64)
--
-- Host: localhost    Database: go
-- ------------------------------------------------------
-- Server version	5.1.61-community

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `species`
--

DROP TABLE IF EXISTS `species`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `species` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `ncbi_taxa_id` int(11) DEFAULT NULL,
  `common_name` varchar(255) DEFAULT NULL,
  `lineage_string` text,
  `genus` varchar(55) DEFAULT NULL,
  `species` varchar(255) DEFAULT NULL,
  `parent_id` int(11) DEFAULT NULL,
  `left_value` int(11) DEFAULT NULL,
  `right_value` int(11) DEFAULT NULL,
  `taxonomic_rank` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `sp0` (`id`),
  UNIQUE KEY `ncbi_taxa_id` (`ncbi_taxa_id`),
  KEY `sp1` (`ncbi_taxa_id`),
  KEY `sp2` (`common_name`),
  KEY `sp3` (`genus`),
  KEY `sp4` (`species`),
  KEY `sp5` (`genus`,`species`),
  KEY `sp6` (`id`,`ncbi_taxa_id`),
  KEY `sp7` (`id`,`ncbi_taxa_id`,`genus`,`species`),
  KEY `sp8` (`parent_id`),
  KEY `sp9` (`left_value`),
  KEY `sp10` (`right_value`),
  KEY `sp11` (`left_value`,`right_value`),
  KEY `sp12` (`id`,`left_value`),
  KEY `sp13` (`genus`,`left_value`,`right_value`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2014-10-01  2:57:56
